#import "@preview/oxdraw:0.1.0": *
#import "@preview/cetz:0.4.2": canvas, draw
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge, shapes



#set document(
  title: "FastDedup: A fast and memory efficient tool for reads deduplication",
  author: "Raphaël Ribes",
  date: auto,
)

#set page(
  paper: "a4",
  margin: (x: 2cm, y: 2cm),
  header: context {
    if counter(page).get().first() > 1 [
      #text(size: 9pt, fill: gray)[FastDedup -- A fast and memory efficient tool for reads deduplication #h(1fr) R.Ribes & C.Mandier]
      #line(length: 100%, stroke: 0.5pt + gray)
    ]
  },
  footer: context {
    let page-num = counter(page).get().first()
    align(center)[#text(size: 9pt)[#page-num]]
  },
)

#set text(size: 11pt, lang: "en")
#set heading(numbering: "1.1")
#set par(justify: true)

// Style des liens
#show link: it => text(fill: blue, it)

// =============================================================================
// Title page
// =============================================================================

#align(center)[
  #v(3cm)
  #text(size: 24pt, weight: "bold")[FastDedup]
  #v(0.5cm)
  #text(size: 18pt, weight: "bold")[A fast and memory efficient tool for reads deduplication]
  #v(3cm)
  #text(size: 12pt)[*Raphaël Ribes*\ Master 2 Bioinformatics - Université de Montpellier]
  
  #text(size: 12pt)[*Céline Mandier*\ Study Engineer - Institut de Science des données de Montpellier]
  #v(1cm)
  #text(size: 11pt, fill: gray)[#datetime.today().display("[day] [month repr:long] [year]")]
  #v(3cm)
]

// =============================================================================
// Abstract
// =============================================================================
#v(1fr)
#rect(fill: luma(245), inset: 1em, width: 100%)[
  *Abstract* -- TODO
]

#pagebreak()
// =============================================================================
// Table of contents
// =============================================================================

#outline(
  title: [Table of contents],
  indent: auto,
  depth: 2,
)

#pagebreak()

// =============================================================================
// Context
// =============================================================================
= -- Context
#v(1em) 
#h(2em) With the rise of New Generation Sequencing (NGS) technologies, the volume of sequencing data has increased exponentially.
With this surge in data, the development of new tools for processing and analyzing sequencing reads was needed.
One of the first and critical step in the analysis pipeline is the deduplication of reads, which involves identifying and removing duplicate sequences that can arise from PCR amplification or sequencing errors.

#v(1em)


#align(center)[
  #figure(
    canvas({
      import draw: *

      // Define the helper function inside the canvas scope
      let draw_read(x, y, h, label_a, label_b) = {
        // Main fragment line
        line((x, y), (x, y + h), stroke: 1.5pt + black)

        // R1 Arrow (Red)
        line((x - 0.2, y), (x - 0.2, y + 0.6), stroke: 1pt + red, mark: (end: ">", fill: red))
        content((x - 0.5, y + 0.3), text(fill: red, size: 8pt)[R1])

        // R2 Arrow (Blue)
        line((x + 0.2, y + h), (x + 0.2, y + h - 0.6), stroke: 1pt + blue, mark: (end: ">", fill: blue))
        content((x + 0.5, y + h - 0.3), text(fill: blue, size: 8pt)[R2])
      }

      // --- Group A ---
      content((0.75, 3.5), text(weight: "bold", size: 1.2em)[A)])
      draw_read(0, 0, 3, "R1", "R2")
      draw_read(1.5, 0, 2, "R1", "R2")
      
      content((0.75, -0.5), text(fill: red, size: 9pt)[$R 1 = R 1$])
      content((0.75, -1.0), text(fill: blue, size: 9pt)[$R 2 != R 2$])

      // --- Group B ---
      content((5.75, 3.5), text(weight: "bold", size: 1.2em)[B)])
      draw_read(5, 0, 3, "R1", "R2")
      draw_read(6.5, 0, 3, "R1", "R2")

      content((5.75, -0.5), text(fill: red, size: 9pt)[$R 1 = R 1$])
      content((5.75, -1.0), text(fill: blue, size: 9pt)[$R 2 = R 2$])
    }),
    caption: [Identifying PCR Duplicates in Paired-End Sequencing.],
  )
]

#v(1em)

#align(center)[
  #block(width: 100%, stroke: .5pt + luma(150), radius: 2pt, inset: 10pt, align(left)[
    #v(0.5em)
    *Solid Black Line:* Represents a single sequenced DNA fragment (height reflects length).
    
    *Red Arrow (#text(fill: red)[R1]):* Forward read (Read 1).
    
    *Blue Arrow (#text(fill: blue)[R2]):* Reverse read (Read 2).
    #v(0.5em)
    *A) Not a duplication:* Shows two distinct fragments where the start position of Read 1 is shared, but different fragment lengths lead to unique mapping positions for Read 2. \
    #v(0.5em)
    *B) PCR Duplication:* Shows two exact fragments where the start position of Read 1 and Read 2 is shared.
  ])
]

We define a fragment $F$, the 5'->3' sequence as $R 1$ with a length of $r_1$ and the 3'->5' sequence as $R 2$ with a length of $r_2$. 
We can say that $R 1$ is the prefix of $F$, $R 2$ reversed is the suffix of $F$ and $r_1 = r_2$.

Therefore, a PCR duplicate can be defined as if $F_1$ and $F_2$ have the same starting and ending position in the genome, with the same length.
In that case, $R 1_1$ and $R 1_2$ are the same, and $R 2_1$ and $R 2_2$ are the same.

This step is crucial because removing these artifacts preventing false-positive during variant calling since PCR polymerases can introduce errors during the cycles of amplification.
If an error occurs early, it is propagated into multiple daughter molecules.
When sequenced, these identical, error-carrying molecules will appear as multiple independent reads supporting a mutation that does not actually exist in the biological sample @Ebbert2016.                                 

Deduplication also allows correcting the allele frequency bias since not all DNA molecules amplify with equal efficiency due to factors like GC content or fragment length.
This uneven amplification can artificially inflate the coverage of certain alleles.
And it is also improving _de novo_ assembly, especially in whole-genome sequencing.
Paired-end reads are mapped to estimate the order and intervening distance between contigs, duplicates can create false-positive connections between contigs or introduce conflicting connections. @Zhang2023.
                             
Before presenting our tool, we will first review the existing tools for reads deduplication and their limitations. 

#pagebreak() 
= -- Existing tools and limitations 
== fastp

`fastp` is a widely adopted, all-in-one FASTQ preprocessor written in C++ that performs quality control, adapter trimming, and deduplication within a single pass @Chen2018 @Chen2023.
Unlike traditional alignment-based tools that require computationally expensive mapping and coordinate sorting, `fastp` operates directly on raw sequences.

To achieve high performance, `fastp` implements a probabilistic deduplication mechanism based on Bloom filters and multiple hash arrays.
For each read (or read pair), the sequence is mapped to multiple integer hash values. 
If all corresponding positions in the hash arrays are already marked as positive, the read is flagged as a duplicate.
This streaming, single-pass architecture avoids reading and writing intermediate files to disk, drastically reducing I/O bottlenecks and resulting in ultra-fast execution speeds.

A significant advantage of `fastp` is its highly configurable memory management.
Users can adjust the `--dup_calc_accuracy` parameter (ranging from level 1 to 6) to scale the hash buffer size from 1 GB up to 24 GB, allowing the tool to run on both standard personal computers and High-performance Computing (HPC) clusters. 

However, this probabilistic approach involves a minor trade-off in absolute precision. `fastp` has a theoretical false-positive (hash collision) rate of approximately 0.01%, which, while negligible for most routine downstream analyses, represents a conceptual difference from exact deterministic methods. 
Furthermore, `fastp`'s deduplication algorithm relies heavily on the exact matching of sequence lengths and coordination regions. 
Consequently, applying certain quality trimming operations (such as sliding-window trimming via `--cut_front` or `--cut_tail`) prior to deduplication can interfere with the duplicate identification process, as trimmed reads may no longer perfectly match their duplicate counterparts.

== FastUniq

`FastUniq` @Xu2012 is a deterministic deduplication tool explicitly designed for paired-end (PE) short reads.
By directly comparing the nucleotide sequences of read pairs, it provides an exact identification of PCR duplicates. 

Algorithmically, `FastUniq` operates through a strict three-step pipeline:
1. loads the entirety of the paired reads into memory
2. sorts them based on their sequences using a merge sort algorithm
3. identifies duplicates by sequentially comparing adjacent read pairs in the sorted list. 

A notable feature of its comparison logic is its ability to handle reads of varying lengths: a shorter sequence is considered identical to a longer one if it perfectly matches its 5' prefix.
If duplicates are found, the pair containing the longest reads is retained.

While `FastUniq` guarantees absolute precision for exact duplicate removal (zero mismatches allowed), its architectural design introduces critical resource bottlenecks. 
Because it relies on importing the entire dataset into RAM simultaneously using a hierarchical storage structure, its memory footprint scales linearly with the input file size.
Benchmarking studies have highlighted this as a severe limitation for modern high-throughput data; for instance, processing a 538 GB Hi-C dataset required up to 1 TB of RAM @Sigorskikh2025. 

Consequently, while highly accurate, `FastUniq` is often impractical for massive datasets without access to HPC clusters. 
Furthermore, its scope is strictly limited to paired-end data, it does not tolerate sequence mismatches, and its deduplication output is known to be sensitive to the order of the input files.

== Clumpify

`Clumpify`, a component of the Java-based BBTools suite @Bushnell2014, introduces a distinct algorithmic paradigm to duplicate removal based on $k$-mer grouping.
Rather than relying on exact sequence hashing or global coordinate sorting, `Clumpify` identifies reads that share common $k$-mers and groups them into localized clusters or "clumps". 

This spatial reorganization provides a unique secondary benefit. 
By physically placing highly similar sequences adjacent to one another in the output file, `Clumpify` dramatically improves standard `gzip` compression efficiency, as the compression algorithm can utilize much closer reference pointers. 

A major functional advantage of `Clumpify` is its high sensitivity to optical duplicates.
By parsing the flowcell tile and pixel coordinates embedded within Illumina FASTQ headers, `Clumpify` (via parameters such as `dupedist`) can explicitly distinguish between amplification-derived PCR duplicates and sensor-derived optical duplicates.
For paired-end data, `Clumpify` predominantly uses the first read (R1) to determine the clump, allowing the paired read (R2) to passively follow its counterpart, which maintains file synchronization while minimizing computational overhead.

To manage memory, `Clumpify` implements a multi-phase strategy (KmerSplit and KmerSort).
If the dataset exceeds available RAM, the tool can divide the data into temporary groups on the disk, sorting them individually before merging.
This approach effectively lowers the time complexity from $O(N log N)$ to approximately $O(N log (N / "groups" ))$, theoretically allowing it to process unlimited data volumes.

However, benchmarking shows that its reliance on the Java Virtual Machine (JVM) can lead to variable stability and high memory overhead on massive datasets.
Sometimes it requires multi-threading just to maintain stability on large inputs like 500+ GB Hi-C datasets @Sigorskikh2025.

== SeqKit

`SeqKit` @Shen2016 @Shen2024 is a comprehensive, cross-platform, and ultra-fast toolkit for FASTA and FASTQ file manipulation.
Written in the Go, it benefits from efficient concurrency management and the ability to produce static binaries with no external dependencies, making it highly accessible and stable across different operating systems.

For PCR deduplication, SeqKit provides the `rmdup` subcommand.
This tool relies on an exact sequence matching approach using a hashing strategy.
To process massive datasets without exhausting system memory, `seqkit rmdup` offers an MD5 mode (invoked via the `-m` or `--md5` flag).
In this mode, the algorithm computes and stores only the MD5 digest of each sequence in the RAM rather than the full nucleotide sequence.
This extreme memory efficiency allows SeqKit to maintain a very low resource footprint, even outperforming several newer tools when handling vast amounts of single-end (SE) data. 

Despite its speed and memory efficiency, `SeqKit rmdup` has a notable architectural limitation: it does not natively support the synchronized deduplication of paired-end reads.
The `rmdup` command evaluates and processes individual files independently.
It will check duplicates inside each file separately, so if $R 1_1=R 1_2$ but $R 2_1 eq.not R 2_2$, $R 1_2$ in will be flagged as a duplicate and removed, while $R 2_2$ will be retained, leading to loss of data and desynchronization of the paired files.
While workarounds exist using other SeqKit subcommands (such as `seqkit pair` to resynchronize files), `SeqKit rmdup` is generally considered the optimal choice strictly for single-end libraries.

== CD-HIT-DUP

`CD-HIT-DUP`, a specialized utility within the broader CD-HIT package @Fu2012 @Huang2010, is designed to identify and remove both identical and nearly-identical duplicates. 
Unlike tools restricted to exact string matching, `CD-HIT-DUP` is particulary robust against sequencing errors. 
It incorporates a greedy incremental clustering algorithm that tolerates a user-defined threshold of mismatches, including both substitutions and insertions/deletions.

`CD-HIT-DUP` relies on a prefix-suffix comparison approach.
The tool first bins reads that share an identical prefix of $k$ nucleotides. 
To accelerate this binning and sorting process, it uses a dual numeric encoding strategy (initially base-5, then base-10), allowing prefixes of up to 27 nucleotides to be compactly represented as 64-bit integers.
Within each bin, sequences are sorted by decreasing length, with the longest sequence acting as the cluster's seed.
The algorithm then evaluates the suffixes of the remaining sequences against this seed, utilizing a $k$-mer word-based index table to statistically predict sequence identity and bypass computationally expensive full pairwise alignments when possible.
For paired-end libraries, `CD-HIT-DUP` conceptually joins the read pairs into a single sequence and checks prefixes at both ends simultaneously to ensure paired-end integrity.

While its error tolerance provides a nuanced, high-sensitivity approach to duplicate removal, `CD-HIT-DUP` introduces a severe computational bottleneck regarding memory consumption.
Because the greedy clustering algorithm must retain representative sequences and cluster bins in RAM, its memory footprint scales aggressively with dataset size and the allowed mismatch threshold.
Benchmarking has shown that processing a relatively modest 50 million single-end read library can consume between 33 GB and 53 GB of RAM.
For massive modern datasets, such as a 538 GB Hi-C dataset, `CD-HIT-DUP` required approximately 1 TB of RAM, often crashing on resource-constrained hardware @Sigorskikh2025. 

Furthermore, practical deployment of `CD-HIT-DUP` in modern automated workflows requires caution.
The tool has known stability issues in containerized environments (e.g., Docker and Singularity), occasionally throwing segmentation faults or I/O assertion errors when inputs are provided via process substitution from compressed files.

== PRINSEQ++

`PRINSEQ++` @Cantu2019 is a highly optimized C++ implementation of the Perl-based `prinseq-lite` program.
It is designed to perform a comprehensive suite of quality control, filtering, trimming, and reformatting tasks for genomic and metagenomic sequence data.
By utilizing POSIX threads (pthreads) for parallelization, it achieves processing speeds more than 10 to 16 times faster than its predecessor while maintaining a relatively low memory footprint. 

For sequence deduplication, `PRINSEQ++` departs from exact string comparison or sorting heuristics by implementing a probabilistic data structure known as a Bloom filter.
During execution, each sequence is transformed using multiple fast, non-cryptographic hash functions, and the corresponding indices in a shared bit-array are flagged.
If a newly evaluated sequence maps entirely to bits that are already set to positive, it is classified as a duplicate.
This probabilistic model is highly advantageous for multi-threaded environments because it enables asynchronous read and write access to the filter with minimal thread blocking. 

Additionally, `PRINSEQ++` natively supports the direct reading and writing of `gzip`-compressed FASTQ and FASTA files without the need for intermediate disk I/O, which significantly reduces hard drive storage requirements.
It also natively supports paired-end reads by processing pairs simultaneously and outputting synchronized R1 and R2 files.

However, the deduplication function (`-derep`) is strictly limited to exact duplicates and does not tolerate sequencing errors.
Furthermore, the probabilistic nature of the Bloom filter introduces a negligible but theoretically non-zero rate of false positives.
Finally, when `PRINSEQ++` is executed with multiple threads, the output sequences are not guaranteed to remain in the same order as the input files, which may require consideration if downstream tools expect strict global index synchronization.

== SAM/BAM ecosystem tools

In addition to tools that operate directly on raw FASTQ sequences, the bioinformatics ecosystem relies heavily on utilities designed for SAM and BAM formats.
PCR deduplication tools like such as `samtools` (via `markdup`) @Li2009 and Picard's `MarkDuplicates` @Picard. 

While these tools are typically used downstream after mapping reads to a reference genome, it is technically possible to convert raw, unaligned FASTQ files directly into unmapped SAM/BAM files to process them through these pipelines.
However, even without the computationally heavy step of genomic alignment, the sheer format conversion from FASTQ to SAM/BAM introduces significant and unnecessary I/O and processing overhead.
Because our primary objective is to evaluate and develop a highly optimized, memory-efficient tool strictly for direct FASTQ-to-FASTQ deduplication, these format-dependent tools fall outside the scope of our comparative analysis.

#page(flipped: true)[
  #figure(
    align(center)[
      #table(
        columns: (auto, auto, auto, auto, auto, 1fr),
        align: (col, row) => if row == 0 { center } else { left },
        stroke: 0.5pt + luma(200),
        fill: (col, row) => if row == 0 { luma(240) } else { none },
        
        [*Tool*], [*Algorithm / Approach*], [*PE Support*], [*Error Tolerance*], [*Memory Usage*], [*Key Limitations*],
        
        [`fastp` @Chen2018 @Chen2023], [Probabilistic (Bloom filter/hashes)], [Yes], [No], [Configurable], [False-positive rate; long runtime],
        
        [`FastUniq` @Xu2012], [Deterministic (Merge sort)], [Yes], [No], [Very High], [Severe RAM bottleneck; sensitive to input file order.],
        
        [`Clumpify` @Bushnell2014], [k-mer grouping], [Yes], [Yes], [Variable], [Variable stability on massive datasets; high memory overhead.],
        
        [`SeqKit rmdup` @Shen2016 @Shen2024], [Exact matching (Hashing/MD5)], [No], [No], [Very Low], [Desynchronizes paired-end files.],
        
        [`CD-HIT-DUP` @Fu2012 @Huang2010], [Greedy incremental clustering], [Yes], [Yes], [Very High], [Stability issues in containerized environments.],
        
        [`PRINSEQ++` @Cantu2019], [Probabilistic (Bloom filter)], [Yes], [No], [Low], [False-positives rate]
      )
    ],
    caption: [Summary of existing tools for reads deduplication and their limitations.],
  )
]

= -- FastDedup: High-Performance FASTX Deduplication
#v(1em)

The primary goal behind the development of `FastDedup` (or `FDedup`) was to create a FASTX PCR deduplication tool that prioritizes maximum speed and memory efficiency. 

To achieve this, the tool relies on `xxh3`, a rapid non-cryptographic hash function, to compute a unique fingerprint for each single-end or paired-end read.
These fingerprints are securely cached in memory using `fxhash`, which provides a low-overhead memory footprint.

FDedup provides users with precise control over the hash collision rate while maintaining high performance. It automatically scales its hashing strategy, seamlessly switching between 64-bit and 128-bit hashes based on the estimated input sequence count and a user-defined collision probability threshold.

Specifically designed for seamless integration into OMICS pipelines, the tool features robust incremental deduplication and auto-recovery mechanisms. 
In the event of an interruption FDedup can safely preload existing hashes from the output file to prevent duplicate processing and smoothly resume operations.
If an uncompressed output file becomes corrupted due to a crash, the tool automatically detects the issue, calculates a fail-safe truncation point, and truncates the file to the last valid sequence before continuing.

Finally, a dry-run mode (`--dry-run` or `-s`) is available to calculate the duplication rate without writing any output files. Because file I/O operations account for the majority of execution time, this mode is exceptionally fast and serves as a highly efficient feature for pipeline planning and resource allocation.

#page(flipped: true)[
    == Architecture
    #figure(
      image("SE.mmd.svg", width: 110%),
      caption: [
        Architectural design of `FastDedup` for Single-End (SE) reads.
      ],
    )
    #v(1fr)
    #figure(
          image("PE.mmd.svg", width: 110%),
          caption: [
            Architectural design of `FastDedup` for Paired-End (PE) reads.
          ],
        )

]

= Benchmarking and performance evaluation

TODO

= Discussion and future perspectives

= Conclusion

= Acknowledgements
#v(1em)

Code architecture and design was developed by Raphaël Ribes, assisted by Gemini 3.1.
Feedback on the algorithm design and ideas for features were provided by Céline Mandier, and the implementation was done by Raphaël Ribes.
The tests suite was mostly developped by Gemini 3.1 and Claude Sonnet 4.6, reviewed by Raphaël Ribes.

This paper was written by Raphaël Ribes, reviewed and edited by Céline Mandier, with the assistance of Gemini 3.1 for grammar, and language correction. 

#v(1fr)
#bibliography("works.bib")