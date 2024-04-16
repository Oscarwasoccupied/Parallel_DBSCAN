# Accelerating Needleman-Wunsch Algorithm for Sequence Alignment

## Abstract

This project aims to design and implement parallel solutions for the genomic sequence alignment problem, leveraging both GPU and multi-core CPU architectures. By adapting the Needleman-Wunsch and Smith-Waterman algorithms for these parallel systems, we expect to achieve significant speedups in aligning genomic sequences. The focus will be on optimizing memory usage, minimizing communication overhead, and effectively balancing the computational load across processing units.

## Background

DNA strands are composed of molecules known as nucleotides, which come in four types: Adenine (A), Thymine (T), Cytosine (C), and Guanine (G). The human genome consists of a long sequence containing millions of these nucleotides in varying orders. Certain segments, which form specific genes, are consistently ordered across individuals. In biological research and applications, comparing two DNA sequences to assess their similarity is of paramount importance. Academically, this comparison can elucidate the phylogenetic relationships of an unidentified DNA sequence or the evolutionary connections between known species. Practically, it helps identify how mutations—such as insertions, deletions, and substitutions—have occurred, pinpointing where the most significant changes along the sequence took place.

Given two sequences of lengths \(N\) and \(M\), a brute-force approach checking every possible nucleotide comparison has a time complexity of \(O(NM)\), which simplifies to \(O(N^2)\) if the sequences are of equal length.

The Needleman-Wunsch algorithm, employing dynamic programming, efficiently aligns two DNA sequences by maximizing the number of matching nucleotides and minimizing the number of gaps introduced. This algorithm, along with its counterpart for local alignment, the Smith-Waterman algorithm, represents the cornerstone of genomic sequence alignment. Both are essential for elucidating genetic relationships, evolutionary histories, and mechanisms of disease progression. The computational intensity of these methods, especially with the exponential increase in available genomic data, underscores the need for optimized solutions.

These algorithms exhibit inherent parallelizability, given the independence in calculating scores for cells within the alignment matrix, despite sequential dependencies (from top, left, and top-left neighbors). This project proposes to dissect the alignment matrix into manageable segments for concurrent processing. This approach leverages data parallelism to tackle the intensive computational demands and task parallelism to efficiently handle dependencies and synchronization, thereby optimizing the alignment process.

### Needleman-Wunsch Algorithm Pseudocode
```plaintext
function NeedlemanWunsch(seqA, seqB, match, mismatch, gap)
  let lenA = length(seqA)
  let lenB = length(seqB)
  initialize matrix[lenA+1][lenB+1]

  // Initialization
  for i from 0 to lenA
    matrix[i][0] = i * gap
  for j from 0 to lenB
    matrix[0][j] = j * gap

  // Matrix filling
  for i from 1 to lenA
    for j from 1 to lenB
      if seqA[i-1] == seqB[j-1]
        score = match
      else
        score = mismatch
      matrix[i][j] = max(matrix[i-1][j-1] + score, // Diagonal
                         matrix[i-1][j] + gap,     // Up
                         matrix[i][j-1] + gap)     // Left

  // Traceback (not shown)

  return matrix
end function
```
This pseudocode outlines the core logic of the Needleman-Wunsch algorithm, focusing on global alignment through dynamic programming. Adaptations and optimizations for parallel execution will be explored to enhance performance on GPU and multi-core CPU platforms.

## The Challenge
A major challenge of this project is the parallelization of dynamic programming algorithms, namely Needleman-Wunsch and Smith-Waterman. These algorithms have inherent sequential dependencies, as the computation of one matrix cell relies on the values of its adjacent cells. Efficiently parallelizing such algorithms is non-trivial, since it involves decomposing the problem in a way that allows for concurrent execution while maintaining the integrity of the final results.

Dynamic programming presents a challenge for parallelism due to:
- **Data Dependencies**: Each step of computation depends on the results of previous steps. Finding parallelism requires identifying independent sub-problems or restructuring the algorithm to minimize inter-task communication.
- **Workload Distribution**: Ensuring an even distribution of the workload among multiple processing units is complex due to the irregular nature of the task dependencies.
- **Memory Access Patterns**: Optimizing memory access to maintain high performance and avoid bottlenecks when different processors access shared data is challenging.

Successfully parallelizing these algorithms on multicore and GPU architectures for genomic sequence alignment could significantly speed up bioinformatics workflows, but doing so requires overcoming substantial parallel computing obstacles. This makes the project a promising opportunity to deepen our understanding of parallel algorithm design and high-performance computing in a scientifically impactful domain.