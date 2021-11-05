import java.io.*;
import java.lang.Math;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.math.BigInteger;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.commons.math3.distribution.BinomialDistribution;

public class MethylBlock {
    // public static final double lambda_gc =  1/35.95704;
    // chr22 Total mC   978946
    // Whole genome Total mCpG/CpG/Genome 11,375,107 / 28,310,100/ 2,835,805,964
    // public static final double lambda_gc =  978946.0/(35200000*0.41); // lambda_gc = 0.06783162
    public static final double lambda_gc =  new BigDecimal(11375107L).divide(new BigDecimal(2835805964L), 10, BigDecimal.ROUND_HALF_UP).doubleValue(); // lambda_gc = 0.009983088
    public static int INTERVAL = 100;
    public static int THRESHOLD = 20;
    public static String chromosome = "";

    public static final double LOG_2 = Math.log(2.0);
    public static final double LOG_10 = Math.log(10.0);

    // numbers greater than 10^MAX_DIGITS_10 or e^MAX_DIGITS_E are considered unsafe ('too big') for floating point operations
    private static final int MAX_DIGITS_10 = 294;
    private static final int MAX_DIGITS_2 = 977; // ~ MAX_DIGITS_10 * LN(10)/LN(2)
    private static final int MAX_DIGITS_E = 677; // ~ MAX_DIGITS_10 * LN(10)

    // chr22   17800433        0.0931677018633541      0.6801722977520527      0.0013460761879122351   true    hetero
    //
    static ReferenceSequenceFile fasta = null;
    public static class Score implements Comparable<Score> {
        int pos;
        double homo_met;
        double hetero_met;
        double homo_no;
        boolean methyl = false;
        char base;
        public Score(int pos_, double homo_met_, double hetero_met_, double homo_no_, boolean methyl_, char b){
            pos = pos_;
            homo_met = homo_met_;
            hetero_met = hetero_met_;
            homo_no = homo_no_;
            methyl = methyl_;
			base = b;
        }
        public int compareTo(Score s){
            return pos - s.pos;
        }
    }

    public static void main(String[] argv){
        PrintWriter out = null;
        PrintWriter wig = null;
        PrintWriter rate_wig = null;
		System.err.println("lambda GC: " + lambda_gc);
        try {
            ArrayList<Score> list = new ArrayList<>();
            File dir = new File(argv[0]);
            chromosome = argv[1];
            if(argv.length > 2){
                INTERVAL = Integer.parseInt(argv[2]);
            }
            if(argv.length > 3){
                fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(argv[3])); // ignore except 1st word
            }
            if(argv.length > 4){
                THRESHOLD = Integer.parseInt(argv[4]);
            }
            if(argv.length > 5){
                if(!argv[5].endsWith(".bed")){
                    System.err.println("The filename for output must have extension '.bed'.");
                    System.exit(-1);
                }
                out = new PrintWriter(new BufferedWriter(new FileWriter(argv[5], true)));
                wig = new PrintWriter(new BufferedWriter(new FileWriter(argv[5].replaceAll(".bed", ".wig"), true)));
                rate_wig = new PrintWriter(new BufferedWriter(new FileWriter(argv[5].replaceAll(".bed", ".rate.wig"), true)));
            }
            wig.println("variableStep chrom=" + chromosome + " span=" + INTERVAL);
            rate_wig.println("variableStep chrom=" + chromosome + " span=" + INTERVAL);
            File[] files = dir.listFiles(
                new FilenameFilter(){
                    public boolean accept(File dir, String name){
                        if(name.startsWith(chromosome + "_") && name.endsWith(".methyl.txt")){
                            return true;
                        }
                        return false;
                    }
                }
            );
            for(File f: files){
                // [0]   [1]      [2] [3]                 [4]                  [5]  [6]         [7]                 [8]                [9]              [10]               [11]
                // chr16 53400046 1.0 0.16666666666666674 0.007936507936507936 true homo_methyl 250,217,184,155,227 18.783606570334165 32.5109375167422 44.901443011937126 c
                // chr16 53400047 0.4666666666666666 0.6083916083916086 0.006993006993006999 true hetero 23,241,248,247,239,247,249,73 30.629296552577113 41.54183738666164 88.86797522551461 g
                //
                // System.out.println(f.getName());
                BufferedReader br = new BufferedReader(new FileReader(f));
                String raw = null;
                while(null != (raw = br.readLine())){
                    String[] line = raw.split("\t");
					if(line.length < 12){
						System.err.println("Bad line: " + raw);
					}
                    // System.out.println(line[0]);
                    list.add(new Score(
                        Integer.parseInt(line[1]), 
                        Double.parseDouble(line[2]), Double.parseDouble(line[3]), Double.parseDouble(line[4]),
                        Boolean.parseBoolean(line[5]), line[11].toUpperCase().charAt(0)));
                }
                br.close();
            }
            java.util.Collections.sort(list);
            for(int i = 0; i<list.size() -1; i++){
                Score s = list.get(i);
                Info info = calcProb(list, i);
                if (info.phred > THRESHOLD && info.gc_rate > 0.1){
                    if(out == null){
                        System.out.println(chromosome + "\t" + s.pos + "\t" + (s.pos + INTERVAL) + "\trate=" + info.gc_rate + ";count=" + info.count + "/" + info.gc_count + ";p="+info.p + ";phred="+ info.phred);
                    }else {
                        out.println(chromosome + "\t" + s.pos + "\t" + (s.pos + INTERVAL) + "\trate=" + info.gc_rate + ";count=" + info.count + "/" + info.gc_count + ";p="+info.p + ";phred="+ info.phred);
                    }
                }
                if(out != null){
                    wig.println(s.pos + "\t" + info.phred);
                    rate_wig.println(s.pos + "\t" + info.gc_rate);
                }
            }
        }catch(Exception e){
            e.printStackTrace();
        }
        finally {
            if(out != null){
                out.close();
                wig.close();
                rate_wig.close();
            }
        }
    }
    public static class Info {
        public int count;
        public int gc_count;
        public double p;
        public double gc_rate;
        public int cpg_count;
        public double phred;
    }
    public static int countCpG(String src){
        int count = 0;
        String seq = src.toUpperCase();
        for(int i = 0; i<src.length()-1; i++){
            if(seq.charAt(i) == 'C' && seq.charAt(i+1) == 'G'){
                count++;
            }
        }
        return count;
    }
    public static Info calcProb(ArrayList<Score> list, int index){
        int count = 0;
		int count_all = 0;
        int start = list.get(index).pos;
        int end = start + INTERVAL;
        for(int i = index; i<list.size()-1; i++){
			Score cur = list.get(i);
			Score next = list.get(i+1);
            if(cur.pos < end && cur.methyl){
                // at least of CpG base was methylated
				if(cur.pos + 1 == next.pos && cur.base =='C' && cur.methyl && next.methyl){ 
					count--;
				}
				count++;
				count_all++;
            }else {
                break;
            }
        }
        Info info = new Info();
        info.count = count;

        String subseq = fasta.getSubsequenceAt(chromosome, start, end).getBaseString();
        info.gc_count = (int)subseq.chars().filter(ch -> (ch == 'C' || ch == 'c' || ch == 'g' || ch == 'G')).count();
        info.cpg_count = countCpG(subseq);
        info.gc_rate = ((double)count)/info.cpg_count;

        // new BinomialDistribution(int trials, double p)
        BinomialDistribution binomial = new BinomialDistribution(info.cpg_count, lambda_gc);

        info.p = 1.0 - binomial.cumulativeProbability(info.count);
        info.phred = -10*Math.log10(info.p);
        if(Double.isInfinite(info.phred)){
            info.phred = 9999.0;
        }
        /*
        double denomitor = logBigInteger(factorialHavingLargeResult(count));
        double numerator = Math.log10(Math.pow(lambda*INTERVAL, count) * Math.exp(-1*lambda*INTERVAL));
        double log_p = -1 * (numerator - denomitor/Math.log(10));

        info.count = count;
        info.log_p = log_p;
        */
        return info;
    }
/*
    public static BigInteger factorialHavingLargeResult(int n) {
        BigInteger result = BigInteger.ONE;
        for (int i = 2; i <= n; i++)
            result = result.multiply(BigInteger.valueOf(i));
        return result;
    }
    public static double logBigInteger(BigInteger val) {
        if (val.signum() < 1)
            return val.signum() < 0 ? Double.NaN : Double.NEGATIVE_INFINITY;
        int blex = val.bitLength() - MAX_DIGITS_2; // any value in 60..1023 works here
        if (blex > 0)
            val = val.shiftRight(blex);
        double res = Math.log(val.doubleValue());
        return blex > 0 ? res + blex * LOG_2 : res;
    }
    */

}
