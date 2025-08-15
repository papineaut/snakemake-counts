rule Heather:
    output:
        "splicingIndex/{sample}/{sample}.splicingIndex.GENCODE_V38.txt"
    params:
        script = "script/junctionCounts_for_splicingIndex.py",
        strandness = config["strandness"]
    input:
        bam="/home/thomaspv/projects/def-kchoquet/files_to_share/iNGN_TGIRTseq_bam/{sample}_Aligned.sortedByCoord.out.bam",
        introns="data/GENCODE_V38_introns.noChr.parsed.bed",

    log:
        "/home/thomaspv/scratch/charles/logs/{sample}.log"
    
    shell:
        """
        module load StdEnv/2020
        module load scipy-stack
        source $HOME/test_env/bin/activate

        python {params.script} {input.introns} {input.bam} {output} {params.strandness}
        """
        