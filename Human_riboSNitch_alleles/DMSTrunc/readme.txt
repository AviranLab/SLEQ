Sequence: Allele_A.fa and Allele_C.fa

DATA: SRR3929636.sam

Candidate structures of allele A and C: Allele_A.txt and Allele_C.txt

PROCESSING:
1. We obtained AlleleReadsAC in DMS-MaPSeq. Therefore we do not need run raw data again.
2. Run DMSTruncation.xcodeproj project, please  initialize your own file path. Two files are saved (AlleleSEQResponse_shortfolds.txt and AlleleSEQDSMatrix_shortfolds.txt).
3. Run NNLS_TRUNC_MaPSeq.R to get the final results. Please initialize your own path.
id: index of structures. 
rhot: corresponding relative abundances