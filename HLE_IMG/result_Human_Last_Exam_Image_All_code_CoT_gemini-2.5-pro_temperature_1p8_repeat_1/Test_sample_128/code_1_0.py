import sys

def solve_chemistry_problem():
    """
    Identifies Compound A from the reaction and provides details.
    """
    # First, check if rdkit is available.
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        rdkit_available = True
    except ImportError:
        rdkit_available = False

    # --- Information based on Chemical Principles ---

    # Define Reactant and Product Information
    aminopyridine = "2-aminopyridine (C5H6N2)"
    phthalaldehyde = "o-phthalaldehyde (C8H6O2)"
    cyanide_source = "HCN (from TMSCN)"  # Simplified for stoichiometry
    
    product_A_name = "2-(pyridin-2-yl)isoindole-1-carbonitrile"
    product_A_formula_text = "C14H9N3"
    byproduct = "Water (H2O)"

    # Print Reaction Information
    print(f"The reaction is a three-component condensation to form Compound A.")
    print(f"Compound A is identified as: {product_A_name}")
    print("\n" + "-"*40)
    print("Overall Balanced Chemical Equation:")
    
    # Stoichiometry is 1:1:1 -> 1:2, as two dehydration steps occur.
    print(f"1 {aminopyridine} + 1 {phthalaldehyde} + 1 {cyanide_source} -> 1 {product_A_formula_text} + 2 {byproduct}")

    print("\n--- Stoichiometric Coefficients ---")
    print(f"Number for 2-aminopyridine: 1")
    print(f"Number for o-phthalaldehyde: 1")
    print(f"Number for Cyanide Source: 1")
    print(f"Number for Compound A ({product_A_formula_text}): 1")
    print(f"Number for Water (byproduct): 2")
    print("-" * 40 + "\n")

    # --- Verification using RDKit (if available) ---
    if rdkit_available:
        print("Verifying properties of Compound A using RDKit...")
        # SMILES string for 2-(pyridin-2-yl)isoindole-1-carbonitrile
        product_smiles = "c1cccc2c1C(C#N)=CN2c1ncccc1"
        mol_A = Chem.MolFromSmiles(product_smiles)
        
        if mol_A:
            calculated_formula = CalcMolFormula(mol_A)
            exact_mass = Descriptors.ExactMolWt(mol_A)
            
            print(f"Molecular Formula: {calculated_formula}")
            print(f"Exact Mass: {exact_mass:.4f}")
        else:
            print("Error: Could not process the molecular structure.", file=sys.stderr)
    else:
        print("Note: RDKit library not found. Cannot programmatically verify properties.")
        print("Please install it for full functionality: pip install rdkit-pypi")

# Run the function
solve_chemistry_problem()
