/*
    GISAID MSA protocol for SARS-CoV-2
    --------------------------------------------------------------------------------------------------------
    This full alignment (msa_0517.fasta) is based on 15543674 submissions to EpiCoV.
    Both duplicate and low-quality sequences (>5% NNNNs) have been removed, using only complete sequences (length >29000 bp). 
    The resulting alignment of 14044655 sequences is created using mafft https://doi.org/10.1093/molbev/mst010 in 3 separate steps.
    
    (1) Each sequence is aligned to the reference hCoV-19/Wuhan/WIV04/2019 (EPI_ISL_402124). Sequences that created dubious insertions 
    of length >=82 nucleotides in the reference sequence or occurred <2 times in the database are not included in the initial alignment. 
    The alignments are created with the command: 
    mafft --addtotop inputBatchOf20kSeqs.fa --thread -1 --6merpair --keeplength --compactmapout refseqWIV04.fa > output.fasta

    (2) Sequences that result in unique insertions in the reference sequence from (1), has occurred >=2 times in the database and contains
    insertions of length <=120, are used as an initial set of sequences for multiple sequence alignment. For this initial set of sequences, 
    we reduce each contiguous stretches of NNNs into a single letter N. This prevents long stretches of NNNNs from causing unnecessarily long 
    insertions in the initial alignment while sequence information carrying true insertions are retained in the alignment. The first sequence 
    in this initial alignment is the reference (EPI_ISL_402124) sequence. The following command is used to align these initial set of sequences 
    to a version of the reference sequence where the first 29 bases are masked by replacing them with NNNs (refseqWIV04Masked29.fa) :
    mafft --addfragments seqsCausingInsertionsInRef.fasta --thread -1 --6merpair refseqWIV04Masked29.fa > seqsCausingInsertionsInRef_aligned.fasta

    (3) The remaining sequences are aligned to the resulting alignment in step (2) with this command:
    mafft --addtotop remainingBatchesOfSequences -thread -1 --6merpair --keeplength --compactmapout reference.fa > msa_0517.fasta

    Note that reference.fa in (3) is composed of 2 sequences. The first is the aligned refseqWIV04Masked29 sequence containing gaps.
    The second is the same aligned reference sequence but the gaps are individually replaced by the letter "N". Mafft version 7.497 which supports 
    the addtotop, compactmapout options is used.

    Acknowledgements
    ================
    We will like to specially thank Kazutaka Katoh (RIMD, Osaka Univ) for formulating this revised multiple alignment strategy and developing new 
    mafft options which significantly improve both the speed and the results of this MSA. We also thank Rob Lanfear for the helpful discussion and suggestions.

    Reference
    =========
    Kazutaka Katoh et. al., NAR 2002 Jul 15;30(14):3059-66. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform
    --------------------------------------------------------------------------------------------------------


*/

/*
    Input
        file of sequences
        request file with reference sequence
        number of inserts to consider
        output file name
    
    Run
        Filter input genomes to get 
            high quality genomes
            remove repeats
        
            Divide the sequence into batches of 20k genomes
            `mafft --addtotop inputBatchOf20kSeqs.fa --thread -1 --6merpair --keeplength --compactmapout reference_sequence.fa > output.fasta`


*/

fn main() {
    println!("Hello, world!");
}
