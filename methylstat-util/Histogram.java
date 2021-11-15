import java.util.zip.GZIPInputStream;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Iterator;

public class Histogram {
	public static class Base {
		public short score;
		public char nuc;
		public Base(short s, char n){
			nuc = n;
			score = s;
		}
	}
	public static Base[] parseInsert(String scores, String bases){
		String[] score_arr = trim(scores).split(", ");
		String[] nuc_arr = trim(bases).split(", ");
		//	for(String s: score_arr){
		//		System.out.println("{" + s + "}");
		//	}
		Base[] buf = new Base[score_arr.length];
		for(int i = 0; i<score_arr.length && score_arr[i].length()>0; i++){
			buf[i] = new Base(Short.parseShort(score_arr[i]),nuc_arr[i].charAt(1));
		}
		return buf;
	}
	public static String trim(String src){
		return src.substring(1, src.length()-1);
	}
	public void run(String path, boolean eval_by_rate, double rate_threshold, boolean methylkit, int modbaseprobs_threshold){
		System.err.println("methylkit mode: " + methylkit);
	    String chromosome = null;
		int start = -1;
		int end =-1;

        File file = new File(path);
        String fname = file.getName();
        String[] buf = fname.split("_");
        chromosome = buf[0];
        String[] buf2 = buf[1].split("-");
        start = Integer.parseInt(buf2[0]);
        end = Integer.parseInt(buf2[1].replaceAll(".txt", ""));

		int length = end - start + 1;
        PrintWriter out = null;
        try {
            out = new PrintWriter(new BufferedWriter(new FileWriter(chromosome + "_" + start + "-" + end + ".methyl.txt")));
        }catch(Exception e){
            System.err.println(e.getMessage());
            System.exit(1);
        }

		ArrayList<LinkedList<Short>> hist = new ArrayList<LinkedList<Short>>();
		ArrayList<LinkedList<Base[]>> insert = new ArrayList<LinkedList<Base[]>>();
		ArrayList<LinkedList<String>> bases = new ArrayList<LinkedList<String>>();
		for(int i = 0; i<length; i++){
			bases.add(new LinkedList<String>());
			hist.add(new LinkedList<Short>());
			insert.add(new LinkedList<Base[]>());
		}
		try {
			BufferedReader br = null;
            if (path.endsWith(".gz")){
                br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
            }else {
                br = new BufferedReader(new FileReader(path));
            }
			String raw = null;
			int[][] table = new int[2][2];
			int lineno = 0;
			while(null != (raw = br.readLine())){
				lineno++;
				// System.err.println("try: " + path);
				// System.out.println(raw);
				// 2629587 @74935df9-c5bb-434d-aee4-4aedb4411a78   0       G       C       -
				String[] line = raw.split("\t");
				// 0: pos, 1: qname, 2: score, 3: genomic, 4: read_base, 5: strand
				int pos = Integer.parseInt(line[0]);
				short score = Short.parseShort(line[2]);
				String genomic_base = line[3];

				hist.get(pos - start +1).add(score);
				bases.get(pos - start + 1).add(genomic_base);
				try {
					insert.get(pos - start + 1).add(parseInsert(line[8], line[9]));
				}catch(ArrayIndexOutOfBoundsException e){
					System.err.println("trancated file: " + path + " at line " + lineno);
					System.exit(0);
				}

			}
			br.close();
			for(int i = 0; i<hist.size(); i++){
				// System.out.print("" + (i + start) + ": ");
				// do exact test
				Iterator<Short> it = hist.get(i).iterator();
				if(hist.get(i).size() == 0){
					continue;
				}
				String base = bases.get(i).get(0);
				int met = 0;
				int nomet = 0;
                double[] bbuf = new double[hist.get(i).size()]; // for Bayes
                short[] data = new short[hist.get(i).size()];
                int j = 0;
				while(it.hasNext()){
					short s = it.next();
					// out.print(s + ", ");
					// System.err.println(path + ": " + s);
					if (s > modbaseprobs_threshold){
						met++;
					}else {
						nomet++;
					}
                    data[j] = s;
                    bbuf[j++] = (s+1)/256.0;
				}
				// out.println("");
        		if(met+nomet == 0 && !methylkit){
        			continue;
				}
                //--Bayes--
                Bayes bayes = new Bayes(bbuf);
                double pMM = bayes.get5mC5mC();
                double pMC = bayes.get5mCC();
                double pCC = bayes.getCC();
                //--Fisher--
				// common 
				table[0][0] = met;
				table[0][1] = nomet;
				// test met homo
				table[1][0] = met+nomet;
				table[1][1] = 0;
				double homo_met = Fisher.exactTest(table);
				table[1][0] = 0;
				table[1][1] = met+nomet;
				double homo_nomet = Fisher.exactTest(table);
				table[1][0] = (met+nomet)/2;
				table[1][1] = (met+nomet)/2;
				double hetero = Fisher.exactTest(table);

				boolean result_met = false;
				String result_type = "no_methyl";
				if(eval_by_rate){
					if(((double)met)/(met + nomet) >= rate_threshold){
						result_met = true;
					}
				}else {
					if (homo_nomet < 0.01){
						result_met = true;
					}
				}
				if(homo_met > hetero && homo_met > homo_nomet){
					result_type = "homo_methyl";
				}else if(hetero > homo_nomet && hetero > homo_met){
					result_type = "hetero";
				}
    			if(!methylkit && result_met == false && result_type == "no_methyl"){
    				continue;
    			}
                String counts = join(data);
				out.print(chromosome + "\t" + (start + i) + "\t" + homo_met + "\t" + hetero + "\t" + homo_nomet +"\t");
				out.println("" + result_met + "\t" + result_type + "\t" + counts + "\t" + pMM + "\t" + pMC + "\t"+ pCC + "\t" + base);
				out.flush();
			}
            out.close();
		}catch(Exception e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
    public static String join(short[] buf){
        StringBuilder s = new StringBuilder();
        s.append(buf[0]);
        for(int i = 1; i<buf.length; i++){
            s.append(",");
            s.append(String.valueOf(buf[i]));
        }
        return s.toString();
    }
}
