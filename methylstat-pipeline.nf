REGION_GEN = '/methylstat-util/region_gen.pl'
METHYLSTAT_BIN = '/methylstat-util/methylstat/debug/methylstat'
CWD = '/nas/test/methyl_ONT_adp_007'
FAI = '/refgenomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa.fai'


process generate_region {

	output:
	stdout into region_list
	"""
	perl ${REGION_GEN} --interval 100000 --fai ${FAI}
	"""
}

process methylstat {
	input:
	tuple val(region), val(name) from region_list.splitText().splitCsv(sep: '\t')
	output:
	file '*.txt' into stat_result

	"""
	$METHYLSTAT_BIN --region $region --index $CWD/fast5index.txt --fasta-index  $FAI --fast5-dir $CWD/workspace --bam $CWD/merged.bam > ${name}.txt
	"""
}

process methylcall {
	publishDir 'called'
	    
	input:
	file statList from stat_result.collect()

	output:
	file '*.methyl.txt' into call_result
	"""
	/methylstat-util/methylcall .
	"""
}

params.CLASSPATH="/methylstat-util:/methylstat-util/lib/gatk-package-4.1.4.1-spark.jar"
params.REFERENCE="/refgenomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa"
params.INTERVAL=500
params.THRESHOLD=20

chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY']


process run_b {
    input:
	val(dummy) from call_result
    val(chr) from chr_list

    script:
    output = "methylblock.interval" + params.THRESHOLD + "." + params.INTERVAL + "." + chr + ".bed"
    """
    echo "java -cp ${params.CLASSPATH} MethylBlock called $chr ${params.INTERVAL} ${params.REFERENCE} ${params.THRESHOLD} ${output}"
    """
}
