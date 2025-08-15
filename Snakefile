configfile: "config.yaml"

include: "rules/Heather.smk"
include: "rules/control.smk"
include: "rules/Stats.smk"

rule all:
    input:
        expand("splicingIndex/{sample}/{sample}.splicingIndex.GENCODE_V38.txt",sample=config['SAMPLES']), # HEATHER 
        expand("results/{comparaison}_rMATS.csv", comparaison=config['COMPARAISONS']), # CONTROL
        #expand("results/test/{comparaison}", comparaison=config['COMPARAISONS']),
        expand("results/{comparaison}/rMATS_Result_P.txt",comparaison=config['COMPARAISONS']), # RMATS-stats1
        expand("results/{comparaison}/rMATS_Result_P_F.txt",comparaison=config['COMPARAISONS']), #RMATS-FDR