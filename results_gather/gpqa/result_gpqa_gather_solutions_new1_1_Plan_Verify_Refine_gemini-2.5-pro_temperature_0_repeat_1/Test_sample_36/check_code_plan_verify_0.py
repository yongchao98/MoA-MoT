from rdkit import Chem
from rdkit.Chem import AllChem

def check_chemistry_answer():
    """
    This function checks the correctness of the answer to the chemistry question.
    It verifies the structure of the final product and counts its 13C-NMR signals.
    """
    # --- Step 1: Verify the structure of the final product (E) ---
    # The reaction sequence is:
    # 1. Propionaldehyde -> 2-ethyl-1,3-dithiolane (A)
    # 2. A + BuLi -> Carbanion (B)
    # 3. B + Bromoethane -> 2,2-diethyl-1,3-dithiolane (C)
    # 4. C -> 3-pentanone (D)
    # 5. D + Wittig reagent from 3-bromopentane -> E
    
    # The ketone (D) is 3-pentanone: (CH3CH2)2C=O. SMILES: CCC(=O)CC
    # The Wittig ylide is formed from 3-bromopentane: (CH3CH2)2CH-Br.
    # The ylide is (CH3CH2)2C=PPh3.
    # The Wittig reaction is: (CH3CH2)2C=O + (CH3CH2)2C=PPh3 -> (CH3CH2)2C=C(CH2CH3)2 + Ph3PO
    # The final product (E) is 3,4-diethylhex-3-ene.
    
    final_product_smiles = "C(CC)(CC)=C(CC)(CC)"
    
    try:
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol is None:
            return f"Error: The SMILES string for the final product '{final_product_smiles}' is invalid."
        
        # --- Step 2: Analyze the symmetry and count 13C-NMR signals ---
        
        # A key point is that this molecule has no E/Z isomers, because the two
        # substituents on each carbon of the double bond are identical (both are ethyl groups).
        # Any answer claiming a mixture of E/Z isomers (like Answer 8) is incorrect.
        
        # Add hydrogens to the molecule for a complete graph representation, which helps
        # in robust symmetry perception.
        mol_with_hs = Chem.AddHs(mol)
        
        # Use canonical atom ranking to find symmetrically equivalent atoms.
        # Atoms with the same rank are considered equivalent.
        ranks = list(Chem.CanonicalRankAtoms(mol_with_hs, breakTies=False))
        
        # Collect the ranks of only the carbon atoms.
        carbon_ranks = set()
        for atom in mol_with_hs.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Atomic number for Carbon
                rank = ranks[atom.GetIdx()]
                carbon_ranks.add(rank)
        
        # The number of unique ranks is the number of 13C-NMR signals.
        calculated_signals = len(carbon_ranks)
        
        # --- Step 3: Compare with the provided answer ---
        # The question's options are A) 11, B) 8, C) 6, D) 3.
        # The provided answer is 'D', which corresponds to 3 signals.
        expected_signals = 3
        
        if calculated_signals == expected_signals:
            return "Correct"
        else:
            return (f"Incorrect. The final product is 3,4-diethylhex-3-ene ({final_product_smiles}). "
                    f"The code calculated {calculated_signals} unique carbon environments, but the correct answer is {expected_signals}. "
                    f"The provided answer is incorrect in its signal count.")

    except ImportError:
        return "Error: RDKit library is not installed. Please install it using 'pip install rdkit'."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)