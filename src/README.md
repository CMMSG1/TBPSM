# TBPSM
Template-based protein structure modeling

To get target PDB file:
    run format : python TBPM.py targetName targetSeq
            or : python TBPM.py .fasta
                    the format of fasta file :
                        > targetName
                        targetSeq
                    ex.
                        >5AOZ 
                        MGSSHHHHHHSSGLVPRGSHMASPVANADVVFDFQNYTAKAGDEVTVDVLVDSKNKPISAMDVKFKVDSPLTIEEIDKES
                        LAFNTTVMTNMAILGANFKSLDDKGEPLVPKDGAAVFTLYVNVPANTPDGTYYVGFNGKNEVHKSNDGSQFTVASKNGAI
                        TVGTPNEEG
    the PDB file is in a folder called targets, and its name is targetName.pdb
                        
                        