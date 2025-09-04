# To run this code, you need to install the rdkit library:
# pip install rdkit
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_mol_features(smiles: str):
    """Analyzes a molecule from a SMILES string and returns its key features."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    # Add hydrogens to count them correctly for proton analysis
    mol_with_hs = Chem.AddHs(mol)
    
    features = {
        "mol": mol,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "has_p_tolyl": mol.HasSubstructMatch(Chem.MolFromSmarts('c1(C)ccc(cc1)-*')),
        "is_epoxide": mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1')),
        "is_ketone": mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][C](=O)[#6]')),
        "methyl_groups": len(mol_with_hs.GetSubstructMatches(Chem.MolFromSmarts('[CH3]'))),
        "aromatic_protons": 0,
        "vinyl_protons": 0,
    }
    
    # Count vinylic and aromatic protons
    for atom in mol_with_hs.GetAtoms():
        if atom.GetAtomicNum() == 6: # Carbon atom
            num_hs = atom.GetTotalNumHs()
            if atom.GetIsAromatic():
                features["aromatic_protons"] += num_hs
            else:
                # Check if it's a vinylic carbon (part of a non-aromatic double bond)
                is_vinylic = any(
                    bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.GetIsAromatic()
                    for bond in atom.GetBonds()
                )
                if is_vinylic:
                    features["vinyl_protons"] += num_hs
                
    return features

def check_correctness():
    """
    Checks the correctness of the provided answer by verifying all constraints and logical steps.
    """
    question_formula = "C11H12O"
    provided_answer_letter = "B"
    
    # --- Define Structures ---
    # Product identified from NMR: 4-(4-methylphenyl)but-3-en-2-one
    product_smiles = "CC(=O)C=CC1=CC=C(C)C=C1"
    
    options = {
        "A": {"name": "2-(1-phenylprop-1-en-2-yl)oxirane", "smiles": "CC(=Cc1ccccc1)C1OC1"},
        "B": {"name": "2-(4-methylstyryl)oxirane", "smiles": "Cc1ccc(C=CC2OC2)cc1"},
        "C": {"name": "2-methyl-3-styryloxirane", "smiles": "CC1OC(C=Cc2ccccc2)C1"},
        "D": {"name": "2-styrylepoxide", "smiles": "c1ccc(C=CC2OC2)cc1"}
    }
    
    errors = []

    # --- Step 1: Verify the Product Structure against NMR data ---
    product_features = get_mol_features(product_smiles)
    # 1H NMR suggests: 2 methyl groups (2x3H), 4 aromatic protons (2x2H), 2 vinyl protons (2x1H)
    # 13C NMR suggests: 1 ketone
    if product_features["methyl_groups"] != 2:
        errors.append(f"Product NMR Mismatch: Deduced product has {product_features['methyl_groups']} methyl groups, but NMR suggests 2.")
    if product_features["aromatic_protons"] != 4:
        errors.append(f"Product NMR Mismatch: Deduced product has {product_features['aromatic_protons']} aromatic protons, but NMR suggests 4.")
    if product_features["vinyl_protons"] != 2:
        errors.append(f"Product NMR Mismatch: Deduced product has {product_features['vinyl_protons']} vinyl protons, but NMR suggests 2.")
    if not product_features["is_ketone"]:
        errors.append("Product NMR Mismatch: Deduced product is not a ketone, but 13C NMR indicates one.")

    # --- Step 2: Check all options against problem constraints ---
    valid_options = []
    for letter, data in options.items():
        features = get_mol_features(data["smiles"])
        if features["formula"] == question_formula:
            valid_options.append(letter)
        else:
            # This option violates the primary constraint.
            print(f"INFO: Option {letter} ({data['name']}) is invalid because its formula is {features['formula']}, not {question_formula}.")

    # --- Step 3: Verify the logic for choosing the correct answer from valid options ---
    if provided_answer_letter not in valid_options:
        errors.append(f"Constraint Violation: The provided answer {provided_answer_letter} is invalid because its molecular formula does not match {question_formula}.")
    
    # The core logic: the starting material must have a p-tolyl group to form the product.
    if not product_features["has_p_tolyl"]:
        errors.append("Reasoning Error: The product structure derived from NMR does not contain a p-tolyl group.")
    
    chosen_option_features = get_mol_features(options[provided_answer_letter]["smiles"])
    if not chosen_option_features["has_p_tolyl"]:
        errors.append(f"Reasoning Error: The chosen answer {provided_answer_letter} does not have the required p-tolyl group to form the product.")

    # Check if any other valid option also has the p-tolyl group (which would make the choice ambiguous)
    for letter in valid_options:
        if letter != provided_answer_letter:
            if get_mol_features(options[letter]["smiles"])["has_p_tolyl"]:
                errors.append(f"Ambiguity: Option {letter} is also a valid isomer with a p-tolyl group, challenging the uniqueness of the answer.")

    # --- Step 4: Verify Reaction Plausibility ---
    if not chosen_option_features["is_epoxide"]:
        errors.append(f"Reaction Mismatch: The chosen starting material {provided_answer_letter} is not an epoxide.")
    if not product_features["is_ketone"]:
        errors.append("Reaction Mismatch: The final product is not a ketone.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. The following issues were found:\n- " + "\n- ".join(errors)

# Run the verification
result = check_correctness()
print(f"\nFinal Verdict: {result}")