rule Stats:
    input:
        "results/{comparaison}_rMATS.csv"
    output:
        #stats="results/test/{comparaison}/rMATS_Result_P.txt",
        stats="results/{comparaison}/rMATS_Result_P.txt"
    conda: "../envs/rmats.yaml"
    shell:
        """
        python script/rMATS_unpaired.py {input} "results/test/{wildcards.comparaison}" 8 0.1 
        """
rule FDR:
    input:
        a="results/{comparaison}/rMATS_Result_P.txt",
        
    output:
        stats="results/{comparaison}/rMATS_Result_P_F.txt",
    params:
        fdr=0.05
    conda: "../envs/rmats.yaml"
    shell:
        """
        python script/FDR.py {input} {output.stats}
        """