rule control:
    input:
        expand("splicingIndex/{sample}/{sample}.splicingIndex.GENCODE_V38.txt", sample=config['SAMPLES'])

    output:
        rmats=expand("results/{comparaison}_rMATS.csv",comparaison=config['COMPARAISONS'])
    params:
        threshold = config['threshold'],
        samples=config['SAMPLES'],
        length=config['read_length'],
        comparaison=config['COMPARAISONS']
    log:
        "logs/control.log"
    script:
        "../script/rMATS.py"
