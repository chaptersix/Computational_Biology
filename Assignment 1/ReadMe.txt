------------------
Usage Instructions
------------------
Development environment: Python 3.7.5, Windows 10


 - alignment_implementation.py implements both algorithms.
    - Input:
            - 3 file paths as arguments: sequence1 sequence2 outputFilePath
    - Output
        - 1 file containing the results of the two  sequences being run through both algorithms

    Note: you must have the python install added to your PATH for this to work.

Usage Examples:
    python alignment_implementation.py pair1\1k4rA_dengue_virus.fasta pair1\5ire_zika_virus.fasta results/pair1.txt
    python alignment_implementation.py pair2\1dzl_HPV.fasta pair2\3mge_HIV.fasta results/pair2.txt
    python alignment_implementation.py pair3\2plv1_polio_virus.fasta pair3\4rhv1_rhino_virus.fasta results/pair3.txt

Precomputed results are in the results folder