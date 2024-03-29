from .models import Snp, target
from django.shortcuts import render
from django.http import HttpResponseRedirect, HttpResponse, JsonResponse
from django import forms
from datetime import datetime
#from django.utils.safestring import mark_safe

from util import get_basenji_sample_info, get_vcf_info, SAD_pipeline_v2, posthoc_p, plot_snp_p, score_box_plot, parse_rsid_v2, valid_rsid #parse_rsid
from util import summary_score
import kipoi
from glob import glob
import os
import json

record_target_pk = {}
record_snp_pk = {}
record_request_rsid = {}
# Create your views here.

class NewTaskForm(forms.Form):
    cell_type = forms.CharField(label="Cell Type", widget=forms.TextInput(attrs={'placeholder': 'pancreas'})) #, help_text='(One cell type, tissue or cell line) <br><br>'
    rsid = forms.CharField(label="rsID", widget=forms.TextInput(attrs={'placeholder': 'rs34584161, rs11708067'})) #, help_text='(Multiple ids are accepted with a separator as commas or spaces) <br><br>'
    #rsid = forms.CharField(label=mark_safe("<br /> rsID or chromosomal position"), widget=forms.TextInput(attrs={'placeholder': 'rs34584161, rs11708067'}), help_text='(Multiple ids are accepted with a separator as commas or spaces)') #

class RsidListForm(forms.Form):
    selected = forms.CharField(widget= forms.TextInput(attrs={'id':'display_value'}))
#     def __init__(self, *args, **kwargs):
#         super(RsidListForm, self).__init__(*args, **kwargs)
#         for rsid in args[0]:
#             self.fields[rsid] = forms.CheckboxInput()

def map_cell(cell_type, login):
    #record_target_pk.append([])
    parse_info = get_basenji_sample_info(cell_type)
    parse_info1 = get_basenji_sample_info(cell_type, by_epi=True)
    existing_tasks = [i.num for i in target.objects.all()]
    for k,v in parse_info1.items():
        for i in v:
            if i not in existing_tasks:
                n = target.objects.create(info=k, num=i, cell=cell_type)
            #print(n.id)
            else:
                n = target.objects.get(num=i)
            record_target_pk[login].append(n.id)

    return list(parse_info.keys())


def index(request):
    login_time = datetime.now().strftime("%M%S%f")

    record_target_pk[login_time] = []
    record_snp_pk[login_time] = []

    if "bs_info" not in request.session:
        request.session["bs_info"] = []
    
    # if this is a POST request we need to process the form data
    if request.method == 'POST' and 'cell_type' in request.body.decode("utf-8") :
        # create a form instance and populate it with data from the request:
        form = NewTaskForm(request.POST)
        form2 = RsidListForm()
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            cell = form.cleaned_data["cell_type"]
            bs_info = map_cell(cell, login_time)
            rsid = form.cleaned_data["rsid"]
            rsids = parse_rsid_v2(rsid) #valid_rsid() #parse_rsid
            print(rsids)
            if len(rsids) ==0:
                # clear bs_info if rsid not available
                bs_info = []
            else:
                #record_snp_pk.append([])
                existing_rsids = [i.rsid for i in Snp.objects.all()]
                vcf_info = get_vcf_info(rsids)
                rsids = []
                for info in vcf_info:
                    annot = None
                    try:
                        chrm, pos, rsid, ref, alt, annot = info
                    except:
                        chrm, pos, rsid, ref, alt = info
                    if annot is not None:
                        rsids.append('('+ annot+') '+rsid)
                    else:
                        rsids.append(rsid)
                    if rsid not in existing_rsids:
                        n = Snp.objects.create(chrm=chrm, pos=pos, rsid=rsid, ref=ref, alt=alt)
                    #record_snp_pk[-1].append(n.id)
                    else:
                        n = Snp.objects.filter(rsid=rsid)
                        if len(n)>1:
                            for ni in n[:len(n)-1]: #cannot use negative indexing here; query set
                                ni.delete()
                        if n[0].chrm=='':
                            n.update(chrm=chrm, pos=pos, rsid=rsid, ref=ref, alt=alt)
                        n = n[0]
                    record_snp_pk[login_time].append(n.id)
            
            #request.session["form"] = {"cell":cell, "rsid":rsid}
            # redirect to a new URL:
            #return HttpResponseRedirect('/thanks/')
            return render(request, 'variants/index.html', {
                "bs_info": bs_info,
                "rsids": rsids,
                "form": form,
                "form2": form2,
                "login": login_time,
                })

    # elif request.method == 'POST':
    #     form = RsidListForm(request.POST)
    #     input(form)
    #     if form.is_valid():
    #         check_rsids = form.cleaned_data["rsids"]
    #         print('-------------------'*3)
    #         print(check_rsids)
    #         print('-------------------'*3)
    #         return render(request, 'variants/predict.html', {
    #             "login": login_time,
    #             "checked_rsids": checked_rsids,
    #             })
        

    # if a GET (or any other method) we'll create a blank form
    else:
        #form = NameForm()
        try: #20201207
            #target.objects.all().delete()
            #Snp.objects.all().delete()
            if len(record_target_pk[login_time])!=0:
                for id in record_target_pk[login_time]:
                    target.objects.get(id=id).delete()
                for id in record_snp_pk[login_time]:
                    Snp.objects.get(id=id).delete()
        except:
            pass
        return render(request, "variants/index.html", {
            "bs_info": [],
            "rsids": [],
            "form": NewTaskForm(),
            "login": login_time,
            })
    
    

def predict(request, login):
    '''
    '''
    #vcf_info = get_vcf_info('data/DM_variant.csv')
    #vcf_info = [rsid]
    
    # 0. check if it is already done before
    #done_list = [item.split('/')[-1].split(' ')[0] for item in glob('variants/static/pred_out/rs*_box.png')]
    #pred_list = [i.rsid for i in Snp.objects.all() if i.rsid not in done_list]
    #all_list = [i.rsid for i in Snp.objects.all()]
    cell_type = target.objects.get(pk=record_target_pk[login][0]).cell
    #all_rsid = [Snp.objects.get(id=i).rsid for i in record_snp_pk[login]]
    if 'selected' in request.POST.keys():
        request_rsids = request.POST['selected'].split(',')
        record_request_rsid[login] = request_rsids
        print('==='*20)
        print(request_rsids)
        print('==='*20)
    all_rsid = [Snp.objects.get(rsid=i).rsid for i in record_request_rsid[login]]
     # 20210524 update login info to design for existed cases of pkl
    #login2 = login+'_'+cell_type+'_'+all_rsid[0]+'_'+str(len(all_rsid))

    found_list = [item for item in glob(f'variants/static/pred_out/{cell_type}_rs*_box.png')]
    found_list_p = [item for item in glob(f'variants/static/pred_out/{cell_type}_rs*_group_p.png')]
    #done_list = []
    #done_list_p = []
    done_dic = {}
    done_dic_p = {}
     
    if len(found_list)>0:
        #done_list = [item for item in found_list if item.split('/')[-1].split('_')[1].split(' ')[0] in all_rsid]
        #done_list_p = [item for item in found_list_p if item.split('/')[-1].split('_')[1].split('>')[0][:-1] in all_rsid]
        for item in found_list:
            found_id = item.split('/')[-1].split('_')[1].split(' ')[0]
            if found_id in all_rsid:
                done_dic[found_id] = item
        for item in found_list_p:
            found_id = item.split('/')[-1].split('_')[1].split('>')[0][:-1]
            if found_id in all_rsid:
                done_dic_p[found_id] = item       
           
    # if all rsids are already existed
    if len(done_dic)==len(all_rsid):
        box_img_path = [done_dic[item].split('/static/')[-1] for item in all_rsid] # done_list] #[glob(f'variants/static/pred_out/{i}*_box.png')[0].split('/static/')[-1] for i in Snp.objects.all()]
        compare_img_path = [done_dic_p[item].split('/static/')[-1] for item in all_rsid]#done_list_p]#[glob(f'variants/static/pred_out/{i}*_group_p.png')[0].split('/static/')[-1] for i in Snp.objects.all()]
        #input(compare_img_path)
        
        # 20210524 update login name to the existing one
        existing_pkl = [filename for filename in glob(f'variants/static/pred_out/*_{cell_type}_{len(all_rsid)}_{all_rsid[0]}_score_1pos.pkl')]
        if len(existing_pkl)>0:
            #login2 = existing_pkl[0].split('/')[-1][:-15]
            dst = f'variants/static/pred_out/{login}_{cell_type}_{len(all_rsid)}_{all_rsid[0]}_score_1pos.pkl'
            os.rename(existing_pkl[0], dst)
            return render (request, "variants/predict.html", {
                "rsid_list": all_rsid,
                "box_img_path": box_img_path, #sorted(
                "compare_img_path": compare_img_path, #sorted(
                "login":login,
                })

    # 1. load vcf info
    #vcf_info = get_vcf_info(all_rsid) #[i.rsid for i in Snp.objects.all()]  
    vcf_info = []
    for rsid in all_rsid:
        info = Snp.objects.get(rsid=rsid)
        vcf_info.append([info.chrm, info.pos, info.rsid, info.ref, info.alt])
    print(vcf_info)

    # 2. load Basenji model for variant effect prediction
    model = kipoi.get_model('Basenji')

    # 3. load Basenji model for single ref->alt snp comparison 
    parse_info = {}
    #for i in target.objects.all():
    for id in record_target_pk[login]:
        i = target.objects.get(id=id)
        k, v = i.info, i.num
        if k not in parse_info.keys():
            parse_info[k] = [v]
        else:
            parse_info[k].append(v)
    
    #for k,v in parse_info.items():
    #    print(f"{k}: {v}")
    out_path = 'variants/static/pred_out'
    #input(f"target id list: {record_target_pk}")
    snp_grp_score = SAD_pipeline_v2(login, vcf_info, parse_info, model, cell_type, out_path=out_path)
    box_img_path = score_box_plot(snp_grp_score, out_path)
    #p_vals, compare_img_path, t_test_p_vals = posthoc_p(snp_grp_score, out_path)
    compare_img_path = posthoc_p(snp_grp_score, out_path)
    #img_path = []
    #for img1, img2 in zip(box_img_path, compare_img_path):
    #    img_path.append(img1)
    #    img_path.append(img2)

    #if len(p_vals)>10:
        #plot_snp_p(p_vals, out_path, login)

    return render(request, "variants/predict.html", {
        #"img_path": img_path,
        "rsid_list": all_rsid,
        "box_img_path": sorted(box_img_path),
        "compare_img_path": sorted(compare_img_path),
        "login":login,
        #"login2":login2
        })

def mutationMap(request, rsid, login):
    '''
    '''
    # 0. check if it is already done before
    cell_type = target.objects.get(pk=record_target_pk[login][0]).cell
    #done_list = [item.split('/')[-1].split(' ')[0] for item in glob('variants/static/pred_out/rs*_multi.png')]
    done_list = glob(f'variants/static/pred_out/{cell_type}_{rsid}_*_multi.png')
    if len(done_list)>0:
        img_paths = [item.split('/static/')[-1] for item in done_list]
        return render (request, f"variants/mutationMap.html", {
            "img_paths": img_paths,
            "login": login
        })

    # 1. load vcf info
    #vcf_info = get_vcf_info([rsid])
    vcf_info = []
    info = Snp.objects.get(rsid=rsid)
    vcf_info.append([info.chrm, info.pos, info.rsid, info.ref, info.alt])

    # 2. load Basenji model for 41 mutagenesis table
    model = kipoi.get_model('Basenji')

    # 3. load Basenji model for single ref->alt snp comparison 
    parse_info = {}
    #for i in target.objects.all():
    for id in record_target_pk[login]:
        i = target.objects.get(id=id)
        k, v = i.info, i.num
        if k not in parse_info.keys():
            parse_info[k] = [v]
        else:
            parse_info[k].append(v)
    
    out_path = 'variants/static/pred_out'
    _ = SAD_pipeline_v2(login, vcf_info, parse_info, model, cell_type, out_path=out_path, type='multi_scan')
    #progress_view(request, vcf_info, parse_info, model, out_path, type='multi_scan')
    #print('snp_grp_score: ', snp_grp_score)
    img_paths = []
    for k in parse_info.keys():
        img_paths.append(f'pred_out/{cell_type}_{rsid}_{k}_multi.png')

    return render (request, f"variants/mutationMap.html", {
        "img_paths": img_paths,
        "login": login
    })

def help(request):
    return render (request, f"variants/help.html")

def about(request):
    return render (request, f"variants/about.html")

def summary(request, login):
    #20210330
    res = summary_score(login)
    plot_snp_p(res, "variants/static/pred_out", login)
    return render (request, f"variants/summary.html", {
        "list": res,
        "login":login,
        "img": f"pred_out/Snps_p_summary_{login}.png"
    })
    #return JsonResponse(res, safe=False)

#def progress_view(request, vcf_info, parse_info, model, out_path, type='multi_scan'):
#    result = SAD_pipeline_v2.delay(vcf_info, parse_info, model, out_path, type)
#    print(result)
#    return render(request, 'display_progress.html', context={'task_id': result.task_id})
