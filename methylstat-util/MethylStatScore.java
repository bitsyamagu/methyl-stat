import java.io.*;
import java.lang.Math;
import java.util.ArrayList;
import java.math.BigInteger;

public class MethylStatScore {
    public static final double lambda =  1/35.95704;
    public static int INTERVAL = 100;
    public static int THRESHOLD = 20;

    public static final double LOG_2 = Math.log(2.0);
    public static final double LOG_10 = Math.log(10.0);

    // numbers greater than 10^MAX_DIGITS_10 or e^MAX_DIGITS_E are considered unsafe ('too big') for floating point operations
    private static final int MAX_DIGITS_10 = 294;
    private static final int MAX_DIGITS_2 = 977; // ~ MAX_DIGITS_10 * LN(10)/LN(2)
    private static final int MAX_DIGITS_E = 677; // ~ MAX_DIGITS_10 * LN(10)

    // chr22   17800433        0.0931677018633541      0.6801722977520527      0.0013460761879122351   true    hetero
    public static class Score implements Comparable<Score> {
        int pos;
        double homo_met;
        double hetero_met;
        double homo_no;
        boolean methyl = false;
        public Score(int pos_, double homo_met_, double hetero_met_, double homo_no_, boolean methyl_){
            pos = pos_;
            homo_met = homo_met_;
            hetero_met = hetero_met_;
            homo_no = homo_no_;
            methyl = methyl_;
        }
        public int compareTo(Score s){
            return pos - s.pos;
        }
    }

    public static void main(String[] argv){
		PrintWriter out = null;
		PrintWriter wig = null;
        try {
            ArrayList<Score> list = new ArrayList<>();
            String chromosome = argv[1];
            File dir = new File(argv[0]);
			if(argv.length > 2){
				INTERVAL = Integer.parseInt(argv[2]);
			}
			if(argv.length > 3){
				THRESHOLD = Integer.parseInt(argv[3]);
			}
			if(argv.length > 4){
				if(!argv[4].endsWith(".bed")){
					System.err.println("The filename for output must have extension '.bed'.");
					System.exit(-1);
				}
				out = new PrintWriter(new BufferedWriter(new FileWriter(argv[4], true)));
				wig = new PrintWriter(new BufferedWriter(new FileWriter(argv[4].replaceAll(".bed", ".wig"), true)));
			}
			wig.println("variableStep chrom=" + chromosome + " span=" + INTERVAL);
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
                // System.out.println(f.getName());
                BufferedReader br = new BufferedReader(new FileReader(f));
                String raw = null;
                while(null != (raw = br.readLine())){
                    String[] line = raw.split("\t");
                    // System.out.println(line[0]);
                    list.add(new Score(
                        Integer.parseInt(line[1]), 
                        Double.parseDouble(line[2]), Double.parseDouble(line[3]), Double.parseDouble(line[4]),
                        Boolean.parseBoolean(line[5])));
                }
                br.close();
            }
            java.util.Collections.sort(list);
            for(int i = 0; i<list.size() -1; i++){
                Score s = list.get(i);
                Info info = calcProb(list, i);
                if (info.log_p > THRESHOLD){
					if(out == null){
                    	System.out.println(chromosome + "\t" + s.pos + "\t" + (s.pos + INTERVAL) + "\tp=" + info.log_p + ";count=" + info.count);
					}else {
                    	out.println(chromosome + "\t" + s.pos + "\t" + (s.pos + INTERVAL) + "\tp=" + info.log_p + ";count=" + info.count);
					}
                }
				if(out != null){
                	wig.println(s.pos + "\t" + info.log_p);
				}
            }
        }catch(Exception e){
            e.printStackTrace();
        }
		finally {
			if(out != null){
				out.close();
				wig.close();
			}
		}
    }
	public static class Info {
		public int count;
		public double log_p;
	}
    public static Info calcProb(ArrayList<Score> list, int index){
        int count = 0;
        int end = list.get(index).pos + INTERVAL;
        for(int i = index; i<list.size(); i++){
            if(list.get(i).pos < end){
                count++;
            }else {
                break;
            }
        }
        double denomitor = logBigInteger(factorialHavingLargeResult(count));
        double numerator = Math.log10(Math.pow(lambda*INTERVAL, count) * Math.exp(-1*lambda*INTERVAL));
        double log_p = -1 * (numerator - denomitor/Math.log(10));

		Info info = new Info();
		info.count = count;
		info.log_p = log_p;
        return info;
    }

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

}
