use File::Basename 'basename', 'dirname';
MAIN:
{
    my @fast5 = `du -ma $ARGV[0]/*fast5|cut -f2`;
    chomp @fast5;
    # print @fast5;
    my $count = 0;
    foreach my $fast5 (@fast5) {
        $count++;
        print STDERR "$count: ", $fast5, "\n";
        my @reads = `h5ls $fast5`;
        my $fname = basename($fast5);
        foreach my $line (@reads) {
            my @line = split(/\s+/, $line);
            print $line[0], "\t", $fname, "\n";
        }
    }
}


