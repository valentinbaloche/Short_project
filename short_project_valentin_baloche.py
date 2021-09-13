""" Assignment of secondary structures from a PDB file.

Usage :
=======
    python short_project_valentin_baloche.py <file>.pdb

"""

__author__= ("Valentin Baloche")
__contact__= ("valentin.baloche@gmail.com")
__date__= "2020-09-14"
__version__= "1.0"


# PACKAGES

import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy import constants as cst
from scipy.spatial import distance


# CLASS

class Residue:
    """ Combination of several atomic and structural characteristics extracted from the PDB file, for a targeted residue.

        Attributes:
            name: name of the residue in multiple letter code (ex: ALA, ARG, etc.)
            position: position of the residue in the primary protein sequence
            n_coord: spatial coordinates (x, y, z) of the main chain's N atom (amine group) of the residue
            h_coord: spatial coordinates (x, y, z) of the main chain's H atom (amine group) of the residue
            c_coord: spatial coordinates (x, y, z) of the main chain's C atom (carbonyl group) of the residue
            o_coord: spatial coordinates (x, y, z) of the main chain's O atom (carbonyl group) of the residue
            turn_status: defines the residue's position the object is likely to create a hydrogen bond with
                         (a '-1' will define the residue that precedes the object, a '+1' the one that follows) 
            secondary_structure: involvement of the residue in a secondary structure (/, helix, beta strand or turn)
    """

    # Instanciation
    def __init__(self, name, position, n_coord, h_coord, c_coord, o_coord):
        """ Initializes the Residue from the keyword arguments.

            Args:
                name: name of the residue in multiple letter code (ex: ALA, ARG, etc.)
                position: position of the residue in the primary protein sequence
                n_coord: spatial coordinates (x, y, z) of the main chain's N atom (amine group) of the residue
                h_coord: spatial coordinates (x, y, z) of the main chain's H atom (amine group) of the residue
                c_coord: spatial coordinates (x, y, z) of the main chain's C atom (carbonyl group) of the residue
                o_coord: spatial coordinates (x, y, z) of the main chain's O atom (carbonyl group) of the residue

                (at the instanciation, turn_status and secondary_structure = '/')
        """

        self.__name = name
        self.__position = position
        self.__n_coord = n_coord
        self.__h_coord = h_coord
        self.__c_coord = c_coord
        self.__o_coord = o_coord
        self.__turn_status = '/'
        self.__secondary_structure = '/'

    # Getter:
    def get_id(self):
        """ Class method which returns the ID of the object.

            Args:
                none

            Returns:
                position and name (format: str(<position>-<name>))
        """

        return str(self.__position) + '-' + self.__name
    
    def get_n_coord(self):
        """ Class method which returns the spatial coordinates of the main chain's N atom (amine group) of the residue.

            Args:
                none

            Returns:
                spatial coordinates (format: tuple(x, y, z))
        """

        return (self.__n_coord[0], self.__n_coord[1], self.__n_coord[2])

    def get_h_coord(self):
        """ Class method which returns the spatial coordinates of the main chain's H atom (amine group) of the residue.

            Args:
                none

            Returns:
                spatial coordinates (format: tuple(x, y, z))
        """

        # For residues which don't have a free amine group (ex: proline)
        if self.__h_coord == 'none':
            return 'none'
        else:
            return (self.__h_coord[0], self.__h_coord[1], self.__h_coord[2])
    
    def get_c_coord(self):
        """ Class method which returns the spatial coordinates of the main chain's C atom (carbonyl group) of the residue.

            Args:
                none

            Returns:
                spatial coordinates (format: tuple(x, y, z))
        """

        return (self.__c_coord[0], self.__c_coord[1], self.__c_coord[2])

    def get_o_coord(self):
        """ Class method which returns the spatial coordinates of the main chain's O atom (carbonyl group) of the residue.

            Args:
                none

            Returns:
                spatial coordinates (format: tuple(x, y, z))
        """

        return (self.__o_coord[0], self.__o_coord[1], self.__o_coord[2])

    def get_turn_status(self):
        """ Class method returning the position of the residue which is likely to create a hydrogen bond with the object Residue (prediction
            based on the calculation of the electrostatic interaction between CO's object and NH's residue).

            Args:
                none

            Returns:
                a positive or negative integer (a '-1' will define the residue that precedes the object, a '+1' the one that follows)
        """

        return self.__turn_status
    
    def get_secondary_structure(self):
        """ Class method which returns the predicted secondary structure the Residue is involved in.

            Args:
                none

            Returns:
                /, helix, beta strand or turn
        """

        return self.__secondary_structure


    # Setter:
    def set_turn_status(self, status):
        """ Class method which allows to change the protected attributes 'turn_status'.

            Args:
                status: positive or negative integer
        """

        self.__turn_status = status

    def set_secondary_structure(self, status):
        """ Class method which allows to change the protected attributes 'secondary_stucture'.

            Args:
                status: /, helix, beta strand or turn
        """

        self.__secondary_structure = status


# FUNCTIONS

def parse_pdb_file(id, pdb_file_with_hydrogens):
    """ Extracts various descriptive elements from a PDB file.

        Args:
            id: the PDB id of the protein sequence
            pdb_file_with_hydrogens: a PDB file after a '-reduce' processing wich allows to determine the hydrogen atoms positions

        Returns:
            protein_sequence: a list of the residues that compose the protein from Nter to Cter
            atoms_name: a list of the atoms that compose the protein from the first N atom to the last atom of the last side chain
            atoms_coord : a list of the coordinates vectors that define the position of each atoms composing the protein
    """

    # Creation of PDBParser object
    p = PDBParser(QUIET = True)

    # Creation of structure object from the PDB file
    structure = p.get_structure(id, pdb_file_with_hydrogens)

    # Elimination of HETATOMS
    model = structure[0]
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                residue_to_remove.append((chain.id, residue.id))
        if len(chain) == 0:
            chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)

    # Collection of protein sequence, atoms name and coordinates
    protein_sequence = []
    atoms_name = []
    atoms_coord = []

    for model in structure:
        for chain in model:
            for residue in chain:
                protein_sequence.append(residue.resname)
                for atom in residue:
                    atoms_name.append(atom.get_name())
                    atoms_coord.append(atom.get_vector())

    return protein_sequence, atoms_name, atoms_coord
              

def harmonizing_vectors(atoms_name, atoms_coord):
    """ Creates a list of lists from atoms_name or atoms_coord vectors to get a final list which is of the same size than the protein sequence.

        Args:
            atoms_name: a list of the atoms that compose the protein from the first N atom to the last atom of the last side chain
            atoms_coord : a list of the coordinates vectors that define the position of each atoms composing the protein

        Returns:
            atoms_name_per_residue: a list of lists of the atoms that compose the protein, the sublist being the atoms composing one residue
            atoms_coord_per_residue: a list of lists of the atoms coordinates, the sublist being the coordinates of the atoms composing one residue
    """

    atoms_name_per_residue = []
    atoms_coord_per_residue = []
    atoms_name_in_residue_i = []
    atoms_coord_in_residue_i = []

    i = 0
    while i < len(atoms_name):

        # We can easily detect the begining of a new residue which starts with a 'N'
        if atoms_name[i] == 'N' and i != 0:
            atoms_name_per_residue.append(atoms_name_in_residue_i)
            atoms_coord_per_residue.append(atoms_coord_in_residue_i)
            atoms_name_in_residue_i = []
            atoms_coord_in_residue_i = []

        atoms_name_in_residue_i.append(atoms_name[i])
        atoms_coord_in_residue_i.append(atoms_coord[i])
        
        if i == len(atoms_name) - 1:
            atoms_name_per_residue.append(atoms_name_in_residue_i)
            atoms_coord_per_residue.append(atoms_coord_in_residue_i)

        i += 1

    return atoms_name_per_residue, atoms_coord_per_residue


def n_coord(i):
    """ Extracts the coordinates of the main chain's N atom for the residue at the position i-1 in the atoms_coord_per_residue list.

        Args:
            i: index + 1

        Returns: vector (x, y, z)   
    """

    if 'N' in atoms_name_per_residue[i-1]:
        return atoms_coord_per_residue[i-1][atoms_name_per_residue[i-1].index('N')]


def h_coord(i):
    """ Extracts the coordinates of the main chain's H atom for the residue at the position i-1 in the atoms_coord_per_residue list.

        Args:
            i: index + 1

        Returns: vector (x, y, z)   
    """

    if 'H' in atoms_name_per_residue[i-1]:
        return atoms_coord_per_residue[i-1][atoms_name_per_residue[i-1].index('H')]
    
    # For residues which don't have a free amine group (ex: proline)
    else :
        return 'none'
    

def c_coord(i):
    """ Extracts the coordinates of the main chain's C atom for the residue at the position i-1 in the atoms_coord_per_residue list.

        Args:
            i: index + 1

        Returns: vector (x, y, z)   
    """

    if 'C' in atoms_name_per_residue[i-1]:
        return atoms_coord_per_residue[i-1][atoms_name_per_residue[i-1].index('C')]
    

def o_coord(i):
    """ Extracts the coordinates of the main chain's O atom for the residue at the position i-1 in the atoms_coord_per_residue list.

        Args:
            i: index + 1

        Returns: vector (x, y, z)   
    """

    if 'O' in atoms_name_per_residue[i-1]:
        return atoms_coord_per_residue[i-1][atoms_name_per_residue[i-1].index('O')]


def protein_synthesis(protein_sequence):
    """ Creates an object Residue for each residue in the protein sequence, and put them sequentially in a list

        Args:
            protein_sequence: a list of the residues that compose the protein from Nter to Cter

        Returns:
            protein: a list of Residue objects  
    """

    i = 1
    protein = []
    for residue in protein_sequence:
        name = Residue(residue, i, n_coord(i), h_coord(i), c_coord(i), o_coord(i))
        protein.append(name)
        i += 1

    return protein


def hbond_energy(CO, NH):
    """Calculate the electrostatic interaction energy between two H-bonding groups.

        Args:
            CO: object Residue = electron acceptor
            NH: object Residue = electron acceptor

        Returns:
            electrostatic energy (in kcal/mole)
    """

    rON = distance.euclidean(CO.get_o_coord(), NH.get_n_coord())
    rCH = distance.euclidean(CO.get_c_coord(), NH.get_h_coord())
    rOH = distance.euclidean(CO.get_o_coord(), NH.get_h_coord())
    rCN = distance.euclidean(CO.get_c_coord(), NH.get_n_coord())

    return 0.084*((1/rON)+(1/rCH)-(1/rOH)-(1/rCN))*332


def hbond_partners(protein):
    """Calculate the electrostatic interaction energy between each Residues to estimate those which
       are most likely to be linked by a hydrogen bond, and modify their 'turn_status' accordingly.

        Args:
            protein: a list of Residue objects
    """

    i = 0
    while i < len(protein):
        energy_vector = []
        j = 0
        while j < len(protein):
            try:
                # The 2 Residues must be spaced of at least 2 residues
                if abs(i-j) > 2:
                    energy_vector.append(hbond_energy(protein[i], protein[j]))
                    j += 1

                else:
                    energy_vector.append(1000)
                    j += 1
            
            # Error avoided: calculation of Hbond for residues which don't have a free amine group (ex: proline)       
            except:
                energy_vector.append(1000)
                j += 1
    
        if min(energy_vector) < - 0.5:
            protein[i].set_turn_status(energy_vector.index(min(energy_vector))-i)

        i += 1


def find_helix(protein):
    """Predicts the position of the helix in the protein sequence, and modify the 'secondary_structure' Residue's attribute accordingly.

        Args:
            protein: a list of Residue objects
    """
    
    i = 0
    while i < len(protein):
        # Are there 2 consecutive residues with an n-turn < 5 ?
        if protein[i].get_turn_status() != '/' and protein[i].get_turn_status() < 5 and protein[i].get_turn_status() > 0:

            # Checking we are not at the end of the sequence so we can test protein[i+1]
            if i != len(protein) - 1:

                if protein[i+1].get_turn_status() != '/' and protein[i+1].get_turn_status() < 5 and protein[i+1].get_turn_status() > 0: 
            
                    # If so, we are in a helix whose size is going to be delimited
                    j = 0
                    while protein[i+j].get_turn_status() != '/' and protein[i+j].get_turn_status() < 5 and protein[i+j].get_turn_status() > 0:
                        protein[i+j].set_secondary_structure('helix')
                        j += 1

                    # n-turn of the last residue to extend the 'helix' status 
                    k = protein[i+j-1].get_turn_status()  
                    while k > 0:
                        protein[i+j-1+k].set_secondary_structure('helix')
                        k -= 1
            
                # Starting back the iteration at the last residue of the helix determined
                    i = i+j-1
        i += 1


def find_beta_strand(protein):
    """Predicts the position of the beta strands in the protein sequence, and modify the 'secondary_structure' Residue's attribute accordingly.

        Args:
            protein: a list of Residue objects
    """

    i = 0
    while i < len(protein):
        # Checking we are not at the end of the sequence so we can test protein[i+1]
        if i != len(protein) - 2 and i != len(protein) - 1:

            if protein[i].get_turn_status() != '/' and protein[i+2].get_turn_status() != '/':
                if protein[i].get_secondary_structure() != 'helix':
                    if protein[i].get_turn_status() - protein[i+2].get_turn_status() == 4 or protein[i].get_turn_status() - protein[i+2].get_turn_status() == 0:
                        protein[i].set_secondary_structure('beta strand')
                        protein[i+1].set_secondary_structure('beta strand')
                        protein[i+2].set_secondary_structure('beta strand')
        i += 1 


def find_turn(protein):
    """Predicts the position of the turns in the protein sequence, and modify the 'secondary_structure' Residue's attribute accordingly.

        Args:
            protein: a list of Residue objects
    """

    for residue in protein:
        if residue.get_turn_status() != '/':
            if residue.get_secondary_structure() != 'helix' and residue.get_secondary_structure() != 'beta strand':
                residue.set_secondary_structure('turn')


def find_dssp_structure(dssp_file):
    """ Extracts the secondary structure prediction from a DSSP file

        Args:
            dssp_file: txt file generated by the DSSP program

        Returns:
            dssp_strucutre: a list of secondary structures assigned for each residue of the protein sequence
    """

    dssp_structure = []
    dssp_residue = []

    with open (dssp_file, 'r') as file:
        for column in file:
            dssp_structure.append(column[16])
            dssp_residue.append(column[13])
    
    # Elimination of the header
    dssp_structure = dssp_structure[28:]
    dssp_residue = dssp_residue[28:]

    # Elimination of the 'gaps' in the sequence  
    i = 0
    while i < len(dssp_structure):
        if dssp_residue[i] == 'X' or dssp_residue[i] == '!':
            del dssp_residue[i]
            del dssp_structure[i]

        i += 1

    # Harmonization with the nomenclature of the program
    i = 0
    while i < len(dssp_structure):

        if dssp_structure[i] == ' ':
            dssp_structure[i] = '/'
        
        if dssp_structure[i] == 'H' or dssp_structure[i] == 'G' or dssp_structure[i] == 'I':
            dssp_structure[i] = 'helix'
        
        if dssp_structure[i] == 'E':
            dssp_structure[i] = 'beta strand'
    
        if dssp_structure[i] == 'T':
            dssp_structure[i] = 'turn'
        
        if dssp_structure[i] != '/' and dssp_structure[i] != 'helix' and dssp_structure[i] != 'beta strand' and dssp_structure[i] != 'turn':
            dssp_structure[i] = '/'

        i += 1

    return dssp_structure

def compare_with_dssp(protein, dssp_structure):
    """ Compares the secondary structures predicted by the program and those predicted by DSSP

        Args:
            protein: a list of Residue objects
            dssp_structure: a list of secondary structures assigned for each residue of the protein sequence

        Prints:
            description of the comparison results
    """

    # Num of structures predicted by the programm
    num_helix = 0
    num_beta_strand = 0
    num_turn = 0

    for residue in protein: 
        if residue.get_secondary_structure() == 'helix':
            num_helix += 1
        
        if residue.get_secondary_structure() == 'beta strand':
            num_beta_strand += 1
        
        if residue.get_secondary_structure() == 'turn':
            num_turn += 1

    num_structures = num_helix + num_beta_strand + num_turn

    # Comparison with structures predicted by DSSP
    correct = 0
    beta_strand_error = 0
    helix_error = 0
    turn_error = 0

    i = 0
    while i < len(protein):
        if protein[i].get_secondary_structure() == 'helix' or protein[i].get_secondary_structure() == 'beta strand' or protein[i].get_secondary_structure() == 'turn':

            if protein[i].get_secondary_structure() == dssp_structure[i]:
                correct+=1
        
            if protein[i].get_secondary_structure() == 'beta strand' and protein[i].get_secondary_structure() != dssp_structure[i]:
                beta_strand_error+=1
    
            if protein[i].get_secondary_structure() == 'helix' and protein[i].get_secondary_structure() != dssp_structure[i]:
                helix_error+=1
    
            if protein[i].get_secondary_structure() == 'turn' and protein[i].get_secondary_structure() != dssp_structure[i]:
                turn_error+=1
    
        i+=1
    
    correct_perc = round(correct/num_structures*100)

    # Structures found by DSSP but not by the program
    miss_helix = 0
    miss_beta_strand = 0
    miss_turn = 0

    i = 0
    while i < len(protein):
        if protein[i].get_secondary_structure() == '/':
            if dssp_structure[i] != '/':
                if dssp_structure[i] == 'helix':
                    miss_helix += 1
                if dssp_structure[i] == 'beta strand':
                    miss_beta_strand += 1
                if dssp_structure[i] == 'turn':
                    miss_turn +=1
        i += 1

    print(f"\n"
    f"====================================== ANALYSIS OF {id}.pdb ====================================== \n"
    f"\n"
    f"{correct_perc}% of the structures predicted in this program are also found in the DSSP prediction \n"
    f"(accuracy for each structure: helix: {num_helix - helix_error}/{num_helix} ; beta strand: {num_beta_strand - beta_strand_error}/{num_beta_strand} ; turn : {num_turn - turn_error}/{num_turn}). \n"
    f"\n"
    f"Comparing the regions predicted as 'unstructured' by this program and the DSSP prediction, \n"
    f"this program has missed {miss_helix} helix, {miss_beta_strand} beta strand(s) and {miss_turn} turn(s).\n"
    f"\n"
    f"For more details, please see the {id}_spvb_results.txt in your working directory.\n"
    f"\n"
    f"==================================================================================================")

def analysis_results(protein, dssp_structure):
    """ Details the secondary structure predicted by the program and DSSP in a dataframe

        Args:
            protein: a list of Residue objects
            dssp_structure: a list of secondary structures assigned for each residue of the protein sequence

        Returns:
            DataFrame: a dataframe containing the id of the residue, the SPVB secondary structure prediction and the DSSP's
    """

    amino_acid = []
    predicted_structure = []

    for residue in protein:
        amino_acid.append(residue.get_id())
        predicted_structure.append(residue.get_secondary_structure())

    return pd.DataFrame({'RESIDUE': amino_acid, '| SPVB-PREDICTION': predicted_structure, '| DSSP-PREDICTION': dssp_structure}, 
                columns = ['RESIDUE', '| SPVB-PREDICTION', '| DSSP-PREDICTION'])


# PROGRAM

# Reduce (to add hydrogens) > PDB file with hydrogens
id = str(sys.argv[1])
id = id.replace('.pdb','')

os.system ('reduce -OH -Quiet ' + sys.argv[1] + ' > ' + id + 'OH.pdb')

# Parse > protein_sequence, atoms_name and atoms_coord
protein_sequence, atoms_name, atoms_coord = parse_pdb_file(id, id+'OH.pdb')

# Harmonize > atoms_name_per_residue, atoms_coord_per_residue
atoms_name_per_residue, atoms_coord_per_residue = harmonizing_vectors(atoms_name, atoms_coord)

# Protein synthesis > protein
protein = protein_synthesis(protein_sequence)

# Identification of hydrogen bonds in the sequence 
hbond_partners(protein)

# Secondary structure prediction
find_helix(protein)
find_beta_strand(protein)
find_turn(protein)

# DSSP prediction > dssp_structure
os.system ('dssp -i ' + id + 'OH.pdb -o ' + id + 'OHdssp.txt')
dssp_structure = find_dssp_structure(id + 'OHdssp.txt')

# Comparison 
compare_with_dssp(protein, dssp_structure)

# Details > DataFrame
results = analysis_results(protein, dssp_structure).to_string()
output_file = open(id+'_spvb_results.txt', 'w')
output_file.write(results)