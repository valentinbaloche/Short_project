# short_project_valentin_baloche
## Subject description
Implementation of the DSSP method for assigning secondary structures (helix and beta strands)

## Bibliography
[Wolfgang Kabsch, Christian Sander. _Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features._ Biopolymers, December 1983; 22(12): 2577-2637.](https://doi.org/10.1002/bip.360221211)

## Clone this repository
```
$ git clone https://github.com/valentinbaloche/short_project.git
```

## Install conda environment
```
$ conda env create -f short_project_valentin_baloche.yml
```

## Running
`$ python short_project_valentin_baloche.py <file>.pdb`

Be carefull to use a <file>.pdb containing a single protein chain.

## Ouputs
The program generates 3 files during the analysis:
- OH.pdb: pdb file with added hydrogens*
- OHdssp.pdb: pdb file with DSSP predicted structures**
- spvb_results.txt: txt file with alignment of both SPVB (Short Project Valentin Baloche) and DSSP predicted structures
  
_*([Reduce](https://github.com/rlabduke/reduce))_ ; _**([DSSP](https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html))_

## Exemple of run: [1bta](https://www.rcsb.org/structure/1BTA)
`$ python short_project_valentin_baloche.py 1bta.pdb`
  
### Shell message:
```
====================================== ANALYSIS OF 1bta.pdb ====================================== 

78% of the structures predicted in this program are also found in the DSSP prediction 
(accuracy for each structure: helix: 37/46 ; beta strand: 14/16 ; turn : 0/3). 

Comparing the regions predicted as 'unstructured' by this program and the DSSP prediction, 
this program has missed 0 helix, 1 beta strand(s) and 11 turn(s).

For more details, please see the 1bta_spvb_results.txt in your working directory.

==================================================================================================

```
### 1bta_spvb_results.txt:
| RESIDUE | SPVB-PREDICTION | DSSP-PREDICTION |
|--------:|----------------:|----------------:|  
| 1-LYS   |    beta strand  |               / |
| 2-LYS   |    beta strand  |     beta strand |
| 3-ALA   |    beta strand  |     beta strand |
| 4-VAL   |    beta strand  |     beta strand |
| 5-ILE   |    beta strand  |     beta strand |
| 6-ASN   |           turn  |     beta strand |
| 7-GLY   |              /  |            turn |
| 8-GLU   |              /  |            turn |
| 9-GLN   |              /  |            turn |
|10-ILE   |              /  |               / |
|11-ARG   |              /  |               / |
|12-SER   |          helix  |               / |
|13-ILE   |          helix  |           helix |
|14-SER   |          helix  |           helix |
|15-ASP   |          helix  |           helix |
|16-LEU   |          helix  |           helix |
|17-HIS   |          helix  |           helix |
|18-GLN   |          helix  |           helix |
|19-THR   |          helix  |           helix |
|20-LEU   |          helix  |           helix |
|21-LYS   |          helix  |           helix |
|22-LYS   |          helix  |           helix |
|23-GLU   |          helix  |           helix |
|24-LEU   |          helix  |           helix |
|25-ALA   |              /  |            turn |
|26-LEU   |              /  |               / |
|27-PRO   |           turn  |               / |
|28-GLU   |              /  |            turn |
|29-TYR   |              /  |            turn |
|30-TYR   |              /  |               / |
|31-GLY   |              /  |               / |
|32-GLU   |              /  |               / |
|33-ASN   |          helix  |               / |
|34-LEU   |          helix  |           helix |
|35-ASP   |          helix  |           helix |
|36-ALA   |          helix  |           helix |
|37-LEU   |          helix  |           helix |
|38-TRP   |          helix  |           helix |
|39-ASP   |          helix  |           helix |
|40-CYS   |          helix  |           helix |
|41-LEU   |          helix  |           helix |
|42-THR   |          helix  |            turn |
|43-GLY   |              /  |            turn |
|44-TRP   |              /  |            turn |
|45-VAL   |              /  |               / |
|46-GLU   |              /  |               / |
|47-TYR   |              /  |               / |
|48-PRO   |              /  |               / |
|49-LEU   |    beta strand  |     beta strand |
|50-VAL   |    beta strand  |     beta strand |
|51-LEU   |    beta strand  |     beta strand |
|52-GLU   |    beta strand  |     beta strand |
|53-TRP   |    beta strand  |     beta strand |
|54-ARG   |    beta strand  |     beta strand |
|55-GLN   |          helix  |               / |
|56-PHE   |          helix  |            turn |
|57-GLU   |          helix  |            turn |
|58-GLN   |          helix  |           helix |
|59-SER   |          helix  |           helix |
|60-LYS   |          helix  |           helix |
|61-GLN   |          helix  |           helix |
|62-LEU   |          helix  |            turn |
|63-THR   |              /  |            turn |
|64-GLU   |              /  |            turn |
|65-ASN   |              /  |               / |
|66-GLY   |          helix  |               / |
|67-ALA   |          helix  |           helix |
|68-GLU   |          helix  |           helix |
|69-SER   |          helix  |           helix |
|70-VAL   |          helix  |           helix |
|71-LEU   |          helix  |           helix |
|72-GLN   |          helix  |           helix |
|73-VAL   |          helix  |           helix |
|74-PHE   |          helix  |           helix |
|75-ARG   |          helix  |           helix |
|76-GLU   |          helix  |           helix |
|77-ALA   |          helix  |           helix |
|78-LYS   |          helix  |           helix |
|79-ALA   |          helix  |           helix |
|80-GLU   |          helix  |            turn |
|81-GLY   |              /  |            turn |
|82-CYS   |              /  |               / |
|83-ASP   |    beta strand  |               / |
|84-ILE   |    beta strand  |     beta strand |
|85-THR   |    beta strand  |     beta strand |
|86-ILE   |    beta strand  |     beta strand |
|87-ILE   |    beta strand  |     beta strand |
|88-LEU   |              /  |     beta strand |
|89-SER   |           turn  |               / |

  
