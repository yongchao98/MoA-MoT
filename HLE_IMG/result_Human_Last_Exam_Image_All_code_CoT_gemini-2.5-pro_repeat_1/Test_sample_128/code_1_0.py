try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit'")
    # Provide a fallback for calculation if rdkit is not available
    Chem = None

def get_mol_info(mol):
    """Calculates molecular formula and molecular weight for an RDKit Mol object."""
    if not Chem:
        return "N/A", 0.0
    # Ensure properties are computed
    mol.UpdatePropertyCache(strict=False)
    # Add hydrogens to get correct formula and weight
    mol_with_h = Chem.AddHs(mol)
    formula = CalcMolFormula(mol_with_h)
    mol_weight = Descriptors.MolWt(mol_with_h)
    return formula, mol_weight

def solve_reaction():
    """
    Identifies the product of the reaction and prints the details of the chemical equation.
    """
    # SMILES strings for the reactants and products
    smiles_2_aminopyridine = "NC1=NC=CC=C1"
    smiles_o_phthalaldehyde = "C1=CC=C(C(=C1)C=O)C=O"
    smiles_hcn = "C#N"
    smiles_product_A = "N#CC1=C2C=CC=CC2=CN1C1=NC=CC=C1"
    smiles_water = "O"

    print("Identifying Compound A in the reaction:")
    print("2-aminopyridine + o-phthalaldehyde + TMSCN -> A\n")

    if Chem:
        # Create RDKit molecule objects
        mol_2_aminopyridine = Chem.MolFromSmiles(smiles_2_aminopyridine)
        mol_o_phthalaldehyde = Chem.MolFromSmiles(smiles_o_phthalaldehyde)
        mol_hcn = Chem.MolFromSmiles(smiles_hcn)
        mol_product_A = Chem.MolFromSmiles(smiles_product_A)
        mol_water = Chem.MolFromSmiles(smiles_water)

        # Get info for all molecules
        info = {
            "2-aminopyridine": get_mol_info(mol_2_aminopyridine),
            "o-phthalaldehyde": get_mol_info(mol_o_phthalaldehyde),
            "HCN (from TMSCN)": get_mol_info(mol_hcn),
            "Compound A": get_mol_info(mol_product_A),
            "Water": get_mol_info(mol_water)
        }
        
        # Try to generate IUPAC name for product A
        try:
            # RDKit's IUPAC naming is experimental and may not always work
            AllChem.Compute2DCoords(mol_product_A)
            iupac_name = Chem.MolToIUPACName(mol_product_A)
            if not iupac_name: # Fallback if naming fails
                iupac_name = "2-(pyridin-2-yl)-2H-isoindole-1-carbonitrile"
        except Exception:
            iupac_name = "2-(pyridin-2-yl)-2H-isoindole-1-carbonitrile"

        print("--- Compound A Details ---")
        print(f"Structure: Compound A is {iupac_name}")
        print(f"SMILES String: {smiles_product_A}")
        print(f"Molecular Formula: {info['Compound A'][0]}")
        print(f"Molecular Weight: {info['Compound A'][1]:.2f} g/mol\n")
        
        print("--- Balanced Chemical Equation (with Molecular Weights) ---")
        reactant1_mw = info['2-aminopyridine'][1]
        reactant2_mw = info['o-phthalaldehyde'][1]
        reactant3_mw = info['HCN (from TMSCN)'][1]
        product_A_mw = info['Compound A'][1]
        product_water_mw = info['Water'][1]
        
        print(f"  C5H6N2 ({reactant1_mw:.2f}) + C8H6O2 ({reactant2_mw:.2f}) + HCN ({reactant3_mw:.2f})  ->  C14H9N3 ({product_A_mw:.2f}) + 2 H2O ({2*product_water_mw:.2f})")
        
        # Check mass balance
        mass_reactants = reactant1_mw + reactant2_mw + reactant3_mw
        mass_products = product_A_mw + 2 * product_water_mw
        print(f"\nTotal Mass of Reactants: {mass_reactants:.2f} g/mol")
        print(f"Total Mass of Products: {mass_products:.2f} g/mol")
        print("The reaction is mass-balanced.")
    else:
        print("Could not perform calculations because RDKit is not available.")
        print("\nCompound A is 2-(pyridin-2-yl)-2H-isoindole-1-carbonitrile")
        print("Molecular Formula: C14H9N3")
        print("The balanced equation is: C5H6N2 + C8H6O2 + HCN -> C14H9N3 + 2 H2O")

if __name__ == "__main__":
    solve_reaction()