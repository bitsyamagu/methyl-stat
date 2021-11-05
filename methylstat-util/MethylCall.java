import java.io.File;
import java.io.FilenameFilter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.CountDownLatch;

public class MethylCall {
	static boolean rate = false;
	static double rate_threshold = 0.2;
	static boolean methylkit = false;
	public static void main(String[] argv){
		try {
			for(int i = 0; i<argv.length-1; i++){
				if(argv[i].equals("--rate")){
					rate = true;
					rate_threshold = Double.parseDouble(argv[i+1]);
				}else if(argv[i].equals("--methylkit")){
					methylkit = true;
				}
			}
			File[] files = new File(argv[argv.length-1]).listFiles(
				new FilenameFilter(){
					public boolean accept(File f, String name){
						if(name.endsWith(".txt")){
							return true;
						}
						else {
							return false;
						}
					}
				}
			);
			ExecutorService pool = Executors.newFixedThreadPool(12);
			CountDownLatch latch = new CountDownLatch(files.length-1);
			for(File file: files){
				pool.submit(new Runnable(){
					public void run(){
						new Histogram().run(file.getPath(), rate, rate_threshold, methylkit);
						synchronized(latch){
						latch.countDown();
						}
						System.err.println("count-down: " + latch.getCount());
					}
				});
			}
			System.err.println("waiting finish");
			latch.await();
			System.err.println("finished");

		//	for(File f: files){
		//		System.out.println(f.getName());
		//	}
		    pool.shutdown();
		

		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
/*
#!/usr/bin/perl -w

use strict;
use File::Basename 'basename';

my $CLASSPATH = "/antares02/analysis/yamagu/methyl/methylstat-java";
my $PATH = $ARGV[0];

my @files = `ls -1 $PATH/chr*txt`;
chomp @files;

my @list = ();
foreach my $file (@files){
    my $name = basename($file);
    my ($chr, $start, $end, $ext) = split(/[_\-\.]/, $name);
    # print $start, ",", $end, "\n";
    my $ref = {'chr' => $chr, 'start' => $start, 'end'=> $end, 'file' => $file};
    push(@list, $ref);
}

my @sorted = sort { $a->{'chr'} cmp $b->{'chr'} || $a->{'start'} <=> $b->{'end'} } @list;

foreach my $s (@sorted){
    # print $s->{'chr'}, ":", $s->{'start'}, "\n";
    my $com = "java -cp $CLASSPATH:$CLASSPATH/gatk-package-4.1.4.1-spark.jar:. Histogram ".$s->{'file'};
    # print $com, "\n";
    print `$com`;
}

*/
