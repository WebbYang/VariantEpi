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
import threading
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
        #pbar.set_description("Processing %s" % rsid)
        #bool_index = ref_vcf['variants/ID']==rsid
        try:
            record = ref_vcf.query('dbsnp.rsid:'+rsid, fields='dbsnp')
            pos = record['hits'][0]['dbsnp']['hg19']['end'] # for snp, either start or end is fine
            chrm, ref, alt = [record['hits'][0]['dbsnp'][item] for item in ['chrom', 'ref', 'alt']]
            res.append([chrm, pos, rsid, ref, alt])           
        except:
            print(f"{rsid} not in refence vcf file.")               
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

def score2logo(df_mat):
    for i in range(df_mat.shape[0]):
        df_mat.loc[i] = rand_transform(df_mat.loc[i].values)
    max_val = df_mat.max().max()
    if max_val>2:
        df_mat = df_mat/max_val*2
    return df_mat

#@shared_task(bind=True)
def SAD_pipeline_v2(extract_vcf, target_id, model, cell_type, out_path, type='single_snp'):
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
        gen_vcf_itvl([item], 'data/subitem')
        if type=='multi_scan':
            rsid_method_map[rsid] = muta_score_v2('data/subitem_itvl.bed', target_id, model)#, [item]
            for k, rsid_map in rsid_method_map[rsid].items():
                fig, (ax1, ax) = plt.subplots(2, 1, figsize=(14,4))
                cbar_ax = fig.add_axes([0.91, 0.1, .03, .4])
                df = pd.DataFrame(rsid_map)
                df = score2logo(df)
                #print(df)
                lm.Logo(df, ax=ax1)

                df = pd.DataFrame(np.array(list(rsid_map.values())),index=rsid_map.keys())                
                df = df/960
                sns.heatmap(df,cmap="RdBu_r",center=0,ax=ax,cbar_ax=cbar_ax)
                ax.set_title(rsid+' '+k)
                ax.set_xlabel('Position around the snp')
                plt.yticks(rotation=0)
                #plt.show()
                fig.tight_layout(rect=[0, 0, .9, 1])
                fig.savefig(f'{out_path}/{rsid}_{k}_multi.png')
                #f.clf()
                plt.close()
            result += i
        else:
            rsid_method_map[rsid] = muta_score_1_pos('data/subitem_itvl.bed', [item], target_id, model, save=False, out_path=out_path)
        #progress_recorder.set_progress(i + 1, total_tasks)
    
    if type=='single_snp':
        with open(f'{out_path}/score_1pos.pkl', 'wb') as h:
            pickle.dump(rsid_method_map, h, protocol=pickle.HIGHEST_PROTOCOL)

        return rsid_method_map
    return result
    

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
        print(list2str([code_to_base(item) for item in batch['inputs'][1, start:start+10]]))
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
    'T':[0,0,0,1]
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
        'CAGE:embryonic': 'expression',
        'CAGE:adult':'expression',
        'CAGE':'expression'
    }

def posthoc_p(score_dic, out_path):
    '''
    parse each rsid:grp:scores to p values 
    '''
    p_dic = {}
    img_list = []
    for rsid in score_dic.keys():
        col_label, col_val = [], []
        for k,v in score_dic[rsid].items():
            col_label.extend(len(v)*[k])
            col_val.extend(v)
        df_test = pd.DataFrame({'score':col_val,'group':col_label})
        rsid = rsid.split(' ')[0]+rsid.split(' ')[1][-3:]
        print(rsid)
        #df_test.to_csv(f'{out_path}/{rsid}_grp_socre.csv') # 20210316
        try:
            posthoc_df = ph.posthoc_conover(df_test, val_col='score', group_col='group', p_adjust = 'holm')
            posthoc_df.to_csv(f'{out_path}/{rsid}_grp_socre.csv') # 20210316
            min_p_val, img_filename = p_val_heatmap(posthoc_df, out_path, rsid)
            img_list.append(img_filename.split('/static/')[-1])
            p_dic[rsid] = min_p_val
        except:
            input('Error: statistical calculation.')
    return p_dic, img_list

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

def plot_snp_p(p_dics, out_path):
    keys, p_list = [], []
    for k, v in p_dics.items():
        if ':' in k:
            k = k.split(' ')[0]
        keys.append(k)
        p_list.append(v)
    plt.bar(np.arange(len(p_list)), p_list)
    plt.xticks(np.arange(len(p_list)),labels=keys,rotation=90)
    plt.ylabel('-log(p) of variant')
    plt.hlines(2,0,len(p_list),linestyles='dashed',colors='cyan')
    plt.title('Most significant difference between chromatin states snp score')
    #plt.show()
    plt.tight_layout()
    plt.savefig(f'{out_path}/Snps_p_summary.png')
    plt.close()

def score_box_plot(snp_grp_score, out_path):
    img_list = []
    for rsid in snp_grp_score.keys():
        plt.title(rsid)
        xticks = list(snp_grp_score[rsid].keys())
        plt.boxplot(x=snp_grp_score[rsid].values(), labels=xticks)
        plt.xticks(rotation=30)
        plt.ylabel('Difference score')
        plt.tight_layout()  
        filename = f'{out_path}/{rsid}_box.png'
        plt.savefig(filename)
        plt.close()        
        img_list.append(filename.split('/static/')[-1])

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
    for rsid in rsids:
        record = ref_vcf.query('dbsnp.rsid:'+rsid, fields='dbsnp')
        if len(record['hits'])>0:
            res.append(rsid)
        else:
            print(rsid+" not found!")
    return res