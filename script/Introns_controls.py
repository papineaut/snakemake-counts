
# In[1]:
print('acces')

import pandas as pd
import numpy as np
import random,sys,re


print('module')

t = snakemake.params.threshold
samples= snakemake.params.samples
out_file = snakemake.output[0]

print( t,samples,out_file)

# Load snoRNA annotations
snodb = pd.read_csv(r"data/snoDB_All_V2.0.tsv",sep=('\t'))

# Load snoRNA host introns
introns_with_snos = pd.read_table(r"data/GENCODE_V38_introns.noChr.parsed.intersect_snodb_v2.bed",
                                  header=None, dtype={0:'str'})

print('files')
introns_with_snos.columns = ["chr_intron","start_intron", "end_intron", "name_intron","score" ,"strand_intron",
                             "chr_sno", "start_sno", "end_sno", "name_sno","ref_sno","strand_sno","overlap"]

# Combine the two, taking only the relevant fields from snodb, and removing duplicates
sno_fields = ["snodb_id", "gene_name","box_type","is_expressed","box_type","conservation_phastcons","host_biotype", "host_gene_name"]
introns_with_snos_nodup_tmp = introns_with_snos[['chr_intron','start_intron','end_intron','strand_intron','chr_sno',
                                                 'start_sno','end_sno','ref_sno']].drop_duplicates().reset_index(drop=True)
introns_with_snos_nodup = introns_with_snos_nodup_tmp.merge(snodb[sno_fields], left_on='ref_sno', right_on='snodb_id', how='left')


# In[ ]:


# Compute length of each intron
introns_with_snos_nodup['intron_length'] = introns_with_snos_nodup['end_intron'] - introns_with_snos_nodup['start_intron']


# In[ ]:


# Make two sets, one of snoRNA host introns, and one of other introns
sno_introns = introns_with_snos_nodup[introns_with_snos_nodup['ref_sno']!='-1'].reset_index(drop=True)
other_introns = introns_with_snos_nodup[introns_with_snos_nodup['ref_sno']=='-1'].reset_index(drop=True)


print('prep files')
s=samples[0]

SI_all=pd.read_csv(f"splicingIndex/{s}/{s}.splicingIndex.GENCODE_V38.txt",sep='\t',dtype={'chr':'str'})
SI_all = SI_all.rename(columns={"5SS_count": f"5SS_count_{s}","3SS_count": f"3SS_count_{s}", "splice_count": f"splice_count_{s}"})

for s in samples[1:]:
    read=pd.read_csv(f"splicingIndex/{s}/{s}.splicingIndex.GENCODE_V38.txt",sep='\t',dtype={'chr':'str'})
    
    SI_all=SI_all.merge(read, on=['gene','intron','chr','start','end','strand'] )
    SI_all = SI_all.rename(columns={"5SS_count": f"5SS_count_{s}","3SS_count": f"3SS_count_{s}", "splice_count": f"splice_count_{s}"})
    #SI_all = SI_all[(SI_all.filter(regex=f"{s}$")>t).any(axis=1)].reset_index(drop=True)  # every samples must have t reads of gene x

print('first for')

# Get splicing index dataframes from control samples in cell line of interest
SI_all['avg_cov']=0
for s in samples:
    SI_all['avg_cov']=SI_all['avg_cov']+(SI_all.filter(regex=f"{s}$")).sum(axis=1)
SI_all['avg_cov']=SI_all['avg_cov']/9
SI_all=SI_all[SI_all['avg_cov']>t]  # the average cov of gene x is above t


D0_D2=pd.DataFrame()
D0_D4=pd.DataFrame()
D2_D4=pd.DataFrame()
def rmats_prep(SI_all,j1,j2):
    
    # GET ID
    D0_D2['ID'] = SI_all[['chr', 'start', 'end', 'strand']].astype(str).agg('_'.join, axis=1)
    
    # GET 
    D0_D2['IJC1']= (SI_all.filter(like=f"SS_count_{j1}_P69").sum(axis=1).astype(str)) + ',' \
     + (SI_all.filter(like=f"SS_count_{j1}_P71").sum(axis=1).astype(str)) +',' + (SI_all.filter(like=f"SS_count_{j1}_P74").sum(axis=1).astype(str))
     
     
    D0_D2['SJC1']= (SI_all.filter(like=f"splice_count_{j1}_P69").sum(axis=1).astype(str)) + ',' \
     + (SI_all.filter(like=f"splice_count_{j1}_P71").sum(axis=1).astype(str)) +',' + (SI_all.filter(like=f"splice_count_{j1}_P74").sum(axis=1).astype(str))
     
     
    D0_D2['IJC2']= (SI_all.filter(like=f"SS_count_{j2}_P69").sum(axis=1).astype(str)) + ',' \
     + (SI_all.filter(like=f"SS_count_{j2}_P71").sum(axis=1).astype(str)) +',' + (SI_all.filter(like=f"SS_count_{j2}_P74").sum(axis=1).astype(str))
     
     
    D0_D2['SJC2']= (SI_all.filter(like=f"splice_count_{j2}_P69").sum(axis=1).astype(str)) + ',' \
     + (SI_all.filter(like=f"splice_count_{j2}_P71").sum(axis=1).astype(str)) +',' + (SI_all.filter(like=f"splice_count_{j2}_P74").sum(axis=1).astype(str))
     
    D0_D2.to_csv(f'{j1}_{j2}_rMATS.csv' ,index=False)
    return(D0_D2)

D0_D2=rmats_prep(SI_all, "D0", "D2")
D0_D4=rmats_prep(SI_all, "D0", "D4")
D2_D4=rmats_prep(SI_all, "D2", "D4")


SI_all = SI_all[['chr','start','end','strand','avg_cov']].drop_duplicates().reset_index(drop=True)


print('2e for')

# In[4]:


# Merge with annotated introns
sno_introns_expr = sno_introns.merge(SI_all, left_on=['chr_intron','start_intron','end_intron','strand_intron'], right_on=['chr','start','end','strand'])
other_introns_expr = other_introns.merge(SI_all, left_on=['chr_intron','start_intron','end_intron','strand_intron'], right_on=['chr','start','end','strand'])
introns_with_snos_nodup_expr=introns_with_snos_nodup.merge(SI_all, left_on=['chr_intron','start_intron','end_intron','strand_intron'], right_on=['chr','start','end','strand'])
genes_with_sno = sno_introns_expr['host_gene_name'].dropna().unique()

# Step 2: Filter all_introns for introns from those genes
introns_from_sno_genes = introns_with_snos_nodup_expr[introns_with_snos_nodup_expr['host_gene_name'].isin(genes_with_sno)]

# In[7]:

print('prep 2')
random.seed(101010)

size_window = 10
expr_window = 10

n_ctrl_introns = 5

def random_intron_selection(ref_id, host_intron_length, host_intron_cov):
    
    # Define size and expression windows (dependent on the intron length)
    # If intron is ~100nt, then the multiplicative factor will be 1 (so +/- 10 nt in size)
    # if intron is ~1000nt, factor will be 10 (so +/- 100 nt in size)
    length_magnitude = host_intron_length / 100
    min_size = host_intron_length - size_window*length_magnitude
    max_size = host_intron_length + size_window*length_magnitude
    
    expr_magnitude = host_intron_cov / 100
    min_expr = int(host_intron_cov - expr_window*expr_magnitude)
    max_expr = int(host_intron_cov + expr_window*expr_magnitude)
    
    # Retrieve all control introns corresponding to these intervals
    possible_introns = other_introns_expr[(other_introns_expr['intron_length']>= min_size) & 
                                         (other_introns_expr['intron_length']<= max_size) &
                                         (other_introns_expr['avg_cov']>= min_expr) & 
                                         (other_introns_expr['avg_cov']<= max_expr)].reset_index(drop=True)
    
    # If host_intron_cov is > 1000 and not enough introns, don't set a maximum because fewer possibilities
    if host_intron_cov >= 1000 and len(possible_introns)<n_ctrl_introns:    
        possible_introns = other_introns_expr[(other_introns_expr['intron_length']>= min_size) & 
                                         (other_introns_expr['intron_length']<= max_size) &
                                         (other_introns_expr['avg_cov']>= 1000)].reset_index(drop=True)
    
    # Check that the number of introns is higher than the number we want to retrieve
    if len(possible_introns) >= 1:
        if len(possible_introns) >= n_ctrl_introns:
            matched_introns_index = random.sample(possible_introns.index.tolist(), k=n_ctrl_introns)
            matched_introns = possible_introns.iloc[matched_introns_index,:]
        elif len(possible_introns) >= 1 and len(possible_introns) < n_ctrl_introns:
            matched_introns = possible_introns
        return(matched_introns)
    else:
        print("Error: no matching introns for ", ref_id)

print('rando')
# In[8]:


# Initiate a list for the control introns
control_introns_list = []

# For each host intron, retrieve up to 5 control introns
for i in range(len(introns_from_sno_genes)):
    ref_id = introns_from_sno_genes.iloc[i]['ref_sno']
    intron_length = introns_from_sno_genes.iloc[i]['intron_length']
    cov = introns_from_sno_genes.iloc[i]['avg_cov']
    
    # Apply function to retrive N control introns
    try:
        matched_introns = random_intron_selection(ref_id, intron_length, cov)
        matched_introns['matched_sno'] = ref_id
        # Add to list
        control_introns_list.append(matched_introns)
    except TypeError:
        continue
        
# Convert list to dataframe
control_introns_df = pd.concat(control_introns_list).drop_duplicates(subset=['chr_intron','start_intron','end_intron']).reset_index(drop=True)


# Write to file for further analysis
control_introns_df[['chr_intron','start_intron','end_intron','strand_intron']].to_csv("intron_ctrl.csv", sep="\t", header=True, index=False)

