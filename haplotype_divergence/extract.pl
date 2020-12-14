open(POSITIONS,"reflist2");
while(<POSITIONS>){
    chomp;
    my ($seqName,$begin,$end) = split(/\s/);
    open(SAMTOOLS,"samtools faidx 1_On_b1v03.fasta $seqName:$begin-$end |");
    while(my $line = <SAMTOOLS>){
    print $line;
    }
    close(SAMTOOLS);
}
close(POSITIONS);
