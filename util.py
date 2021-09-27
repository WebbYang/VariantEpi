from functools import total_ordering
import kipoi
import allel
from tqdm import tqdm 
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import os
import scikit_posthocs as ph
import myvariant
import logomaker as lm
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from glob import glob
from matplotlib import cm
from scipy.stats import ttest_ind
#import threading
#import json
#from celery import shared_task
#from celery_progress.backend import ProgressRecorder


def get_basenji_sample_info(keyword, by_epi=False): #, ret_identifier=False
    '''
    keyword[string]: cell type (or epi)
    ret_identifier: whether to return sample identifier or not 
    output[dictionary]: basenji_cell_line of k: description/epi, v: number 
    '''
    basenji_cell_line = {} # k: description/epi, v: number 
    #if ret_identifier:
    #    identifier = []
    '''with open('data/sample_wigs.txt','r') as h:
        for i,line in enumerate(h):
            desc = line.split('\t')[-1]
            if keyword in desc:
                print('%d: %s'%(i,desc))
                if ret_identifier:
                    identifier.append(line.split('\t')[0])
                if by_epi:
                    info = desc.split(':')
                    if info[0]!='HISTONE':
                        cell_type = info[1].split(' ')[0]
                        epi = info[0]
                    else:
                        cell_type = info[1].split(' ')[1]
                        epi = info[1].split(' ')[0]
                    if epi not in basenji_cell_line:
                        basenji_cell_line[epi] = [i]
                    else:
                        basenji_cell_line[epi].append(i)
                else:
                    basenji_cell_line[desc] = [i]
    if ret_identifier:
        return basenji_cell_line, identifier'''

    docs = []
    with open('data/sample_wigs.txt','r') as h:
        for line in h:
            desc = line.split('\t')[-1]
            docs.append(desc)
        
    vectorizer = TfidfVectorizer()
    tf_idf = vectorizer.fit_transform(docs)
    qtf_idf = vectorizer.transform([keyword])
    res = cosine_similarity(tf_idf, qtf_idf)
    sims = res.ravel() 
    for i in np.arange(len(res))[sims>0.2]:
        print("{}\t{:.3f}".format(docs[i].rstrip(), res[i][0]))
        if by_epi:
            info = docs[i].split(':')
            if info[0]!='HISTONE':
                cell_type = info[1].split(' ')[0]
                epi = info[0]
            else:
                cell_type = info[1].split(' ')[1]
                epi = info[1].split(' ')[0]
            if epi not in basenji_cell_line:
                basenji_cell_line[epi] = [i]
            else:
                basenji_cell_line[epi].append(i)
        else:
            basenji_cell_line[docs[i].rstrip()] = [i]  

    return basenji_cell_line

def get_vcf_info(file_path, save=False):
    '''
    input type 1. read snp rsid from vcf file column one 
    input type 2. read pos, alt ... to do 
    return vcf info: list of [chr, location, rsid, ref, alt]
    '''
    if type(file_path)==list:
        snps = file_path
    elif file_path.split('.')[-1]=='pkl':
        with open(file_path, 'rb') as f:
            res = pickle.load(f)
        return res
    elif file_path.split('.')[-1]=='tsv':
        snp_df = pd.read_csv(file_path, sep='\t')
        snps = snp_df['SNP'].tolist()
    elif file_path.split('.')[-1]=='csv':
        snp_df = pd.read_csv(file_path)
        snps = snp_df['variant'].tolist()
    else:
        print("Please provide a csv or tsv file.")
        return -1
    #ref_vcf = allel.read_vcf(ref_vcf)
    ref_vcf = myvariant.MyVariantInfo()

    res = []

    pbar = tqdm(snps) # debug [:1]
    for rsid in pbar:
        POSITION_FLAG = False
        #pbar.set_description("Processing %s" % rsid)
        #bool_index = ref_vcf['variants/ID']==rsid
        if ':' in rsid: # 18:5956883 -> chr18:5956883-5956884
            if 'chr' not in rsid:
                rsid = f'chr{rsid}'
            if '-' not in rsid:
                start = int(rsid.split(':')[1])
                rsid = f'{rsid}-{start}' # 'start+1' will be zero-based 
            POSITION_FLAG = True
        try:
            #record = ref_vcf.query('dbsnp.rsid:'+rsid, fields='dbsnp')
            record = ref_vcf.query(rsid, scopes='dbsnp.rsid')
            #pos = record['hits'][0]['dbsnp']['hg19']['end'] 
            # for snp, either start or end is fine, but may encouter indel            
            for i, hit in enumerate(record['hits']):
                if POSITION_FLAG:
                    if hit['dbsnp']['hg19']['end']==start:
                        pos = start
                        idx = i
                        break
                else:
                    pos = hit['dbsnp']['hg19']['start']
                    if hit['dbsnp']['hg19']['end']==pos:
                        idx = i
                        break

            chrm, ref, alt, rsid = [record['hits'][idx]['dbsnp'][item] for item in ['chrom', 'ref', 'alt', 'rsid']]
            if alt=='':
                alt = 'N'
            #if len(ref)>1:
            #    ref = #start
            annot = None
            # for item in record['hits']:
            #     if 'cadd' in item.keys():
            #         annot = item['cadd']['consequence']
            #         if type(annot)==list:
            #             annot = annot[-1]
            #         break
            found = record['hits'][0]['snpeff']['ann']
            if type(found)==list:
                for ann in found:
                    if 'effect' in ann.keys():
                        annot = ann['effect']
                        break
            else:
                if 'effect' in found.keys():
                    annot = found['effect']

            if annot is not None:
                res.append([chrm, pos, rsid, ref, alt, annot])
            else:
                res.append([chrm, pos, rsid, ref, alt])           
        except:
            print(f"{rsid} not in the dbsnp or it's a structural variant.")               
    if save:
        prefix = file_path.split('/')[-1].split('.')[0]
        print(f'Save {prefix}.pkl at data/')
        with open(f'data/{prefix}.pkl', 'wb') as f:
            pickle.dump(res, f)
    return res

##========= newly add 20201130 ==============

def gen_vcf_itvl(vcf_list, filename):
    '''
    Generate interval file format for Basenji prediction.
    chr1	2496649	2496659
    '''
    h = open(filename+'_itvl.bed','w')

    for item in vcf_list:
        str_test = ''
        try:
            chromsome, start, end = 'chr'+item[0], int(item[1])-4, int(item[1])+4
            str_test+='%s\t%d\t%d'%(chromsome,start,end)
            str_test+='\n'
            h.write(str_test)
        except:
            print(item, len(item))
    h.close()  

def rand_transform(row):
    '''
    Mutation map conservation score calculation.
    To be replaced by a formal definition from reference...
    '''
    zero_idx = (row>=-0.0001) & (row<=0.0001)
    new_row = np.zeros(4)
    new_row[zero_idx] = np.log2(-sum(row[row<0])+1.5)
    return new_row

def ref_transform(row):
    zero_idx = row==0
    new_row = np.zeros(4)
    new_row[zero_idx] = 1
    return new_row

def pn_transform(row):
    zero_idx = row==0
    new_row = row.copy()
    #new_row[zero_idx] = 2
    return new_row

def score2logo(df_mat, trans_type = pn_transform):
    for i in range(df_mat.shape[0]):
        df_mat.loc[i] = trans_type(df_mat.loc[i].values)

    return df_mat#/np.abs(df_mat.max()-df_mat.min())*2

#@shared_task(bind=True)
def SAD_pipeline_v2(login, extract_vcf, target_id, model, cell_type, out_path, type='single_snp'):
    '''
    input: 
        extract_vcf:[chr, location, rsid, ref, alt], 
        target_id: (k, v) pair from basenji output
        type: single_snp or multi_scan
    output: sad score around each snp
    '''
    rsid_method_map = {}
    # build uniq rows by rsid
    uniq_rsid_vcf = {}
    for item in extract_vcf:
        if type=='single_snp':
            k = f"{item[2]} chr{item[0]}:{item[1]}{item[3]}>{item[4]}" #  {item[0]}:{item[1]} | {item[0]}:{item[1]}{item[3]}>{item[4]}
        else:
            k = item[2]
        if k not in uniq_rsid_vcf.keys():
            uniq_rsid_vcf[k] = item
            # if there's particular snp to alter
            if k.split(' ')[0]=='rs13266634':
                uniq_rsid_vcf[k][4] = 'A'

    # Calculate score map
    total_tasks = len(list(uniq_rsid_vcf))
    #progress_recorder = ProgressRecorder(self)
    result = 0
    for i,(rsid,item) in enumerate(uniq_rsid_vcf.items()):
        rsid = cell_type+'_'+rsid # to add cell type info to the heading of filename
        print(item)
        #rsid = item[2]
        # Generate inteveral file
        gen_vcf_itvl([item], f'data/subitem')
        if type=='multi_scan':
            rsid_method_map[rsid] = muta_score_v2(f'data/subitem_itvl.bed', target_id, model)#, [item]
            for k, rsid_map in rsid_method_map[rsid].items():
                #fig, (ax1, ax) = plt.subplots(2, 1, figsize=(14,4))
                fig, (ax2, ax1, ax) = plt.subplots(3, 1, figsize=(15,5), gridspec_kw={'height_ratios': [0.4,2.5,2]})
                cbar_ax = fig.add_axes([0.91, 0.1, .03, .3])
                df = pd.DataFrame(rsid_map)
                ref_logo = lm.Logo(score2logo(df, trans_type = ref_transform), ax=ax2)
                ref_logo.ax.set_xticks([])
                ref_logo.ax.set_yticks([])
                ax2.axis("off")
                #ref_logo.ax.text(40, 0.6, 'Reference')
                ref_logo.ax.text(40.5, 0.5, 'Reference', fontsize=14)

                df = score2logo(pd.DataFrame(rsid_map))
                cons_logo = lm.Logo(df, ax=ax1)
                cons_logo.ax.set_xticks([])
                #cons_logo.ax.set_ylim([0, 2])
                df = pd.DataFrame(np.array(list(rsid_map.values())),index=rsid_map.keys())                
                #df = df/960
                sns.heatmap(df,cmap="RdBu_r",center=0,ax=ax,cbar_ax=cbar_ax)
                ax.set_title(rsid+' '+k)
                ax.set_xlabel('Position around the snp')
                plt.yticks(rotation=0)
                #plt.show()
                #fig.tight_layout(rect=[0, 0, .9, 1])
                fig.savefig(f'{out_path}/{rsid}_{k}_multi.png')
                #f.clf()
                plt.close()
            result += i
        else:
            rsid_method_map[rsid] = muta_score_1_pos(f'data/subitem_itvl.bed', [item], target_id, model, save=False, out_path=out_path)
        #progress_recorder.set_progress(i + 1, total_tasks)
    
    if type=='single_snp':
        all_rsids = [item[2] for item in extract_vcf]
        with open(f'{out_path}/{login}_{cell_type}_{len(all_rsids)}_{all_rsids[0]}_score_1pos.pkl', 'wb') as h:
            pickle.dump(rsid_method_map, h, protocol=pickle.HIGHEST_PROTOCOL)

        return rsid_method_map
    return result

#============20210330=============
def find_max_func(dic_data, df_pvals):
    '''
    input: key: 'Open chromatin', 'Enhancer' ..., value: scores
    output: max value key and value
    '''   
    max_func = ''
    max_score = 0
    ttest_p = df_pvals['ttest_p']
    df_pvals = df_pvals.iloc[:,:-1]
    func_list = [item for item in df_pvals.index]
    for k in func_list:
        score = np.mean([abs(item) for item in dic_data[k]])
        if score > abs(max_score):
            max_score = np.mean(dic_data[k]) #dic_data[k][scores.index(max(scores))] 
            max_func = k
            max_p = df_pvals[k].max()
    
    # 20210914 replace
    # max_func = ''
    # max_score = 0
    # ttest_p = df_pvals['ttest_p']
    # df_pvals = df_pvals.iloc[:,:-1]
    # df_v = df_pvals.max()
    # max_p = df_v.max()
    # func_list = [item for item in df_v[df_v==df_v.max()].index]
    # # 20210626: to add the whole vector of p value to compare the function consistency
    # #sum_p = {}
    # for k in func_list:
    #     scores = [abs(item) for item in dic_data[k]]
    #     score = dic_data[k][scores.index(max(scores))]
    #     #sum_p[k] = sum(df_pvals[k])
    #     if abs(score)>abs(max_score):
    #         max_score = score
    #         max_func = k

    # 20210626: check if thae max function contains the max sum p
    # if max_p!=0:
    #     sum_p_reverse = {v:k for k,v in sum_p.items()}
    #     max_sum_p = max(sum_p.values())
    #     if sum_p[max_func] < max_sum_p:
    #         print('The max score is not the max function!')
    #         print(f'original max func: {max_func}')
    #         print(f'original max score: {max_score}')
    #         max_func = sum_p_reverse[max_sum_p]
    #         scores = [abs(item) for item in dic_data[max_func]]
    #         max_score = dic_data[max_func][scores.index(max(scores))]       
    #         max_p = df_v[max_func]
    #         print(f'new max func: {max_func}')
    #         print(f'new max score: {max_score}')
    #         #input('checkpoint')
    # else:
    #     max_func = 'N/A'

    '''for k, scores in dic_data.items():
        for score in scores:
            if abs(score)>abs(max_score):
                max_score = score
                max_func = k
    max_p = df_pvals[max_func].max()'''
    max_ttest_p = ttest_p[max_func]
    #if max_p < 2:
    if max_ttest_p < 2:
        max_p = "NA"
        max_func = "NA"
    else:
        #max_p = '{:.1f}'.format(max_p)
        mask = np.array(func_list)!=max_func
        max_p = df_pvals[max_func][mask].min()
        second_func = np.array(func_list)[mask][df_pvals[max_func][mask]==max_p][0]
        max_p = str(second_func) if max_p<1 else 'NA'
        #print(f'Second function: {second_func}')
    return max_func, '{:.3f}'.format(max_score), max_p, max_ttest_p

def summary_score(login):
    filename = glob(f"variants/static/pred_out/{login}*_score_1pos.pkl")[0]
    with open(filename, 'rb') as h:
        rsid_method_map = pickle.load(h)
    summary_scores = []
    for k, snp_scores in rsid_method_map.items():
        # 20210402
        for file in glob("variants/static/pred_out/"+k.split(' ')[0]+'*_group_p.csv'):
            df_pvals = pd.read_csv(file, index_col=0)

        max_func, max_score, second_func, ttest_p = find_max_func(snp_scores, df_pvals)
        cell, rsid, pos = re.split('_| ', k)
        summary_scores.append((cell, rsid, pos ,max_func, max_score, ttest_p, second_func)) 
    
    return sorted(summary_scores, key=lambda x: float(x[-3]), reverse=True) #abs(float(x[-2]))
#============20210330=============

def muta_score_v2(itvl_file, target_map, model): #, vcf_info
    '''
    itvl_file: bed file generated depend on vcf position
    vcf_info: [chr, location, rsid, ref, alt]
    target_map: (k, v) pair from Basenji ouotput
    '''
    test_kwargs = {'intervals_file': itvl_file,
                   'fasta_file': 'data/hg19.fa'}
    # Get the dataloader and instantiate it
    dl_test = model.default_dataloader(**test_kwargs)
    # get a batch iterator
    it = dl_test.batch_iter(batch_size=1)
    start = int(131072/2)-21

    #for i,item in enumerate(vcf_info):
        #print(i)
    batch = next(it)
    # initial batch, later will use again and again for mutagenesis
    init_batch = batch['inputs']
    
    # Generate reference score
    batch['inputs'] = np.append(batch['inputs'],batch['inputs'], axis=0) # make it two 
    ref_score = model.predict_on_batch(batch['inputs'])
    
    # Build score table
    method = {k:{'A':[], 'T':[], 'C':[], 'G':[]} for k in target_map.keys()}
    #method = {k:{'A':np.zeros(41), 'T':np.zeros(41), 'C':np.zeros(41), 'G':np.zeros(41)} for k in target_map.keys()}
    
    for pos_idx in range(41):
        batch['inputs'] = np.append(init_batch, init_batch, axis=0) # make it four
        batch['inputs'] = np.append(batch['inputs'], batch['inputs'], axis=0)
        for i,base in enumerate(['A','T','C','G']):
            batch['inputs'][i,start+pos_idx] = base_to_code[base]
        
        alt_score1 = model.predict_on_batch(batch['inputs'][:2])
        alt_score2 = model.predict_on_batch(batch['inputs'][2:])
        
        for k,v in target_map.items():
            ref_i = ref_score[0,:,v]
            method[k]['A'].append(np.sum(alt_score1[0,:,v]-ref_i))
            method[k]['T'].append(np.sum(alt_score1[1,:,v]-ref_i))            
            method[k]['C'].append(np.sum(alt_score2[0,:,v]-ref_i))
            method[k]['G'].append(np.sum(alt_score2[1,:,v]-ref_i))

    return method

def muta_score_1_pos(itvl_file, vcf_info, target_map, model, out_path, save=True):
    '''
    itvl_file: bed file generated depend on vcf position
    vcf_info: [chr, location, rsid, ref, alt]
    target_map: (k, v) pair from Basenji ouotput
    ''' 
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    test_kwargs = {'intervals_file': itvl_file,
                   'fasta_file': 'data/hg19.fa'}
    # Get the dataloader and instantiate it
    dl_test = model.default_dataloader(**test_kwargs)
    # get a batch iterator
    it = dl_test.batch_iter(batch_size=1)
    start = int(131072/2)-5
    # store score by rsid:grp_key:[...] 
    score_table = {}
    global group_map

    for i,item in enumerate(vcf_info):
        #print(i)
        batch = next(it)
        rsid, ref, alt = item[2], item[3], item[4]
        # initial batch, later will use again and again for mutagenesis
        init_batch = batch['inputs']
        
        # Generate reference score
        batch['inputs'] = np.append(batch['inputs'],batch['inputs'], axis=0) # make it two   
        batch['inputs'][1,start+4] = base_to_code[alt]
        print(list2str([code_to_base(item) for item in batch['inputs'][0, start:start+10]]))
        #print(list2str([code_to_base(item) for item in batch['inputs'][1, start:start+10]]))
        result = model.predict_on_batch(batch['inputs'])

        # store target score in rsid:grp:[scores]
        score_table = {}

        # find max differential expression from the cell-type related target
        max_score, max_i = 0, 0
        for epi,item in target_map.items():
            try:
                grp = group_map[epi]
            except:
                print(f"Missing experiment type: {epi}")
                continue
            if grp not in score_table.keys():
                score_table[grp] = []
            for i in item:
                score = np.sum(result[1,:,i]-result[0,:,i])
                score_table[grp].append(score)
                #print("%s: %.3f" %(map_sample_wig_desc(i),score))
                if np.abs(score)>np.abs(max_score):
                    max_score, max_i = score, i
        print('%s: %.1f'%(map_sample_wig_desc(max_i), max_score))
        if save:
            file_name = rsid+'_'+alt
            with open(f"{out_path}/exp_profile_"+file_name+".pkl", "wb") as fp: 
                pickle.dump(result, fp, protocol=pickle.HIGHEST_PROTOCOL)

    return score_table

def map_sample_wig_desc(idx):
    with open('data/sample_wigs.txt','r') as h:
        all_rows = h.readlines()
    return all_rows[idx].split('\t')[-1].rstrip()

def code_to_base(code):
    if code[0]:
        return 'A'
    if code[1]:
        return 'C'
    if code[2]:
        return 'G'
    if code[3]:
        return 'T'

def list2str(list_data):
    a = ''
    for item in list_data:
        a+=item
    return a

base_to_code = {
    'A':[1,0,0,0],
    'C':[0,1,0,0],
    'G':[0,0,1,0],
    'T':[0,0,0,1],
    'N':[0,0,0,0]
}

group_map = {'DNASE': 'Open chromatin',
        'H3K9me3': 'Repressed',
        'H3K4me1': 'Enhancer',
        'H3K36me3': 'Active',
        'H3K79me2': 'Active',
        'H3K27ac': 'Enhancer',
        'H3K27me3': 'Repressed',
        'H3K4me3': 'Promoter',
        'H3K9ac': 'Promoter',
        'CAGE:embryonic': 'Expression',
        'CAGE:adult':'Expression',
        'CAGE':'Expression'
    }

def posthoc_p(score_dic, out_path):
    '''
    parse each rsid:grp:scores to p values 
    '''
    p_dic = {}
    t_test_p_dic = {}
    img_list = []
    for rsid in score_dic.keys():
        col_label, col_val = [], []
        for k,v in score_dic[rsid].items():
            col_label.extend(len(v)*[k])
            col_val.extend(v)           
        df_test = pd.DataFrame({'score':col_val,'group':col_label})
        # abs_score = df_test['score'].apply(abs)
        # max_grp = df_test['group'][abs_score==max(abs_score)].iloc[0]
        # print('========')
        # print(max_grp)
        df_test.to_csv('test_34584161.csv')
        # print('========')
        rsid = rsid.split(' ')[0]+rsid.split(' ')[1][-3:]
        # 20210912 add t test p value
        ttest_p_val = []
        for k in ['Active','Enhancer','Expression','Open chromatin','Promoter','Repressed']:
            case = df_test['score'][df_test['group']==k] #max_grp
            control = df_test['score'][df_test['group']!=k] #max_grp
            ttest_p_val.append(ttest_ind(case, control).pvalue)
        #t_test_p_dic[rsid] = np.log(p_val)
        #df_test.to_csv(f'{out_path}/{rsid}_grp_socre.csv') # 20210316
        #try:
        posthoc_df = ph.posthoc_conover(df_test, val_col='score', group_col='group', p_adjust = 'holm')
        min_p_val, img_filename = p_val_heatmap(posthoc_df, out_path, rsid)
        posthoc_df['ttest_p'] = ttest_p_val
        log_df = -np.log10(posthoc_df)
        log_df.to_csv(f'{out_path}/{rsid}_group_p.csv', float_format='%.1f') # 20210316
        
        img_list.append(img_filename.split('/static/')[-1])
            #p_dic[rsid] = min_p_val
        #except:
            #print(f'Error: statistical calculation in {rsid}.')
    #return p_dic, img_list, t_test_p_dic
    return img_list

def p_val_heatmap(posthoc_df, out_path, rsid=None):
    #mask = np.zeros(posthoc_df.shape, dtype=bool)
    #mask[np.triu_indices(len(mask))] = True
    mask = np.eye(posthoc_df.shape[0], dtype=bool)
    data = -np.log10(posthoc_df)
    maxd = data.max().max()
    sns.heatmap(data, vmin = 0, vmax = maxd if maxd>=2 else 2, center = 0, cmap = 'coolwarm', annot = True, mask = mask)
    plt.xticks(rotation = 30)
    if rsid is not None:
        plt.title('Functional change significance pairwise test of \n'+rsid)
    else:
        plt.title('Functional change significance pairwise test')
    #plt.show()
    plt.tight_layout()
    #rsid = rsid.split(' ')[0]+rsid.split(' ')[1][-3:]
    img_filename = f'{out_path}/{rsid}_group_p.png'
    #save_inner = 'variants/templates/variants/'+img_filename
    plt.savefig(img_filename)
    plt.close()
    return data.max().max(), img_filename

def plot_snp_p(res, out_path, login):
    '''
    cell, rsid, pos ,max_func, max_score, max_p, ttest_p
    '''
    rsids = np.array([item[1] for item in res])
    pvals = np.array([float(item[-2]) if item[-1]!='No significance' else 0 for item in res])
    grps = [item[3] for item in res]
    grps = np.array([item if item!='expression' else 'Expression' for item in grps])
    scores = np.array([float(item[-3]) for item in res])

    rsids = rsids[pvals>=2]
    p_list = pvals[pvals>=2]
    grps = grps[pvals>=2]
    scores = scores[pvals>=2]

    c_map = {'Open chromatin':0,'Promoter':1,'Enhancer':2,'Active':3, 'Repressed':4, 'Expression':5, 'N/A':6}
    r_cmap = {v:k for k,v in c_map.items()}
    colors = cm.get_cmap('hsv')
    x = np.arange(len(rsids))
    y = np.array(scores)
    c = np.array([c_map[g] for g in grps])
    #s = np.array([p*5 for p in pvals])       
    #p_list = np.array(pvals)
        
    fig, ax = plt.subplots()
    
    #ax.bar(np.arange(len(p_list)), p_list)
    ax.bar(np.arange(len(p_list)), y)
    plt.xticks(np.arange(len(p_list)),labels=rsids,rotation=90)
    ax.set_ylabel('gain/loss of function')
    ax.bar(np.arange(len(p_list))[y>0], y[y>0])
    plt.title('Most Significant Difference score of each Variant')
    plt.tight_layout() 
    
    ax2 = ax.twinx()
    for i in range(6):
        if len(y[c==i])!=0:
            ax2.scatter(x[c==i], p_list[c==i], c=colors(c[c==i]/6), label=r_cmap[i]) #, s=s[c==i]
    ax2.legend() 
    ax2.hlines(2,-1,len(p_list),linestyles='dashed',colors='k')#colors='cyan'
    ax2.set_ylabel('-log(p) of variant')
    ax2.set_ylim(0, max(p_list)+0.5)
    plt.xlim(-1,len(p_list))
    #plt.show()
    plt.tight_layout()
    plt.savefig(f'{out_path}/Snps_p_summary_{login}.png')
    plt.close()

def old_plot_snp_p(res, out_path, login):
    rsids = [item[1] for item in res]
    pvals = [float(item[-2]) if item[-2]!='No significance' else 0 for item in res]
    grps = [item[3] for item in res]
    grps = [item if item!='expression' else 'Expression' for item in grps]
    
    scores = [float(item[4]) for item in res]

    c_map = {'Open chromatin':0,'Promoter':1,'Enhancer':2,'Active':3, 'Repressed':4, 'Expression':5, 'N/A':6}
    r_cmap = {v:k for k,v in c_map.items()}
    colors = cm.get_cmap('hsv')
    x = np.arange(len(rsids))
    y = np.array(scores)
    c = np.array([c_map[g] for g in grps])
    s = np.array([p*5 for p in pvals])
        
    p_list = np.array(pvals)
        
    fig, ax = plt.subplots()
    
    ax.bar(np.arange(len(p_list)), p_list)
    plt.xticks(np.arange(len(p_list)),labels=rsids,rotation=90)
    ax.set_ylabel('-log(p) of variant')
    plt.hlines(2,-1,len(p_list),linestyles='dashed',colors='k')#colors='cyan'
    plt.title('Most Significant Difference score of each Variant')
    plt.tight_layout()
    
    ax.bar(np.arange(len(p_list))[y>0], p_list[y>0])
    ax2 = ax.twinx()
    for i in range(6):
        if len(y[c==i])!=0:
            ax2.scatter(x[c==i], y[c==i], c=colors(c[c==i]/6), label=r_cmap[i]) #, s=s[c==i]
    ax2.legend() 
    ax2.set_ylabel('gain/loss of function')
    plt.xlim(-1,len(p_list))
    #plt.show()
    plt.tight_layout()
    plt.savefig(f'{out_path}/Snps_p_summary_{login}.png')
    plt.close()

def score_box_plot(snp_grp_score, out_path):
    img_list = []
    for rsid in snp_grp_score.keys():
        plt.title(rsid)
        xticks = list(snp_grp_score[rsid].keys())
        xticks = [item if item!='expression' else 'Expression' for item in xticks]
        plt.boxplot(x=snp_grp_score[rsid].values(), labels=xticks)
        plt.xticks(rotation=30)
        plt.ylabel('Differential score')
        plt.tight_layout()  
        filename = f'{out_path}/{rsid}_box.png'
        plt.savefig(filename)
        plt.close()        
        img_list.append(filename.split('/static/')[-1])
        filename = f'{out_path}/{rsid}_box.csv'
        outdic = {}
        for k,v in snp_grp_score[rsid].items():
            outdic[k] = ["{:.3f}".format(item) for item in [min(v), np.quantile(v,0.25), np.median(v), np.quantile(v,0.75), max(v)]]
        
        df_out = pd.DataFrame(outdic, index=['Min', '25%', '50%', '75%', 'Max'])
        df_out.to_csv(filename)

    return img_list

def parse_rsid(line):
    if "," in line:
        data = line.split(",")
        return [item.strip() for item in data]
    return [line]

def parse_rsid_v2(line):
    if "," in line:
        data = line.split(",")
        return [item.strip() for item in data]
    if " " in line:
        data = line.split(" ")
        return [item for item in data]
    return [line]

def valid_rsid(rsids):
    res = []
    ref_vcf = myvariant.MyVariantInfo()

    # deal with non rsid query, ex: genes
    for rsid in rsids:
        if rsid[:2]!='rs':
            rsids.remove(rsid)
            record = ref_vcf.query(rsid, scopes='dbsnp.rsid')
            for item in record['hits']:
                if 'dbsnp' in item.keys():
                    if 'rsid' in item['dbsnp']:
                        res.append(item['dbsnp']['rsid'])


    for rsid in rsids:
        #record = ref_vcf.query('dbsnp.rsid:'+rsid, fields='dbsnp')
        record = ref_vcf.query(rsid, scopes='dbsnp.rsid')
        #if len(record['hits'])>0:
        try:
            res.append(record['hits'][0]['dbsnp']['rsid'])
        except:
            print(rsid+" not found!")
    return res