# To run this code, you need to have the RDKit library installed.
# You can install it using pip: pip install rdkit-pypi

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def get_molecule_info(name, smiles):
    """Creates an RDKit molecule object and returns its properties."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        return {"name": name, "smiles": smiles, "formula": formula, "mw": mw}
    else:
        return {"name": name, "error": "Invalid SMILES string"}

def main():
    """Main function to identify and describe Compound A."""
    # Define the molecules involved in the overall reaction
    reactant1 = get_molecule_info("3-hydroxy-pyridine-2-carbaldehyde", "O=Cc1c(O)cccn1")
    reactant2 = get_molecule_info("Aniline", "Nc1ccccc1")
    reactant3_name = "Sodium Cyanide" # NaCN, ionic, no SMILES representation for the salt as a whole
    
    # Based on the reaction mechanism, this is the final product A
    # The reaction is an imine formation followed by a Strecker-type addition of cyanide.
    # Structure: 2-(cyano(phenylamino)methyl)pyridin-3-ol
    compound_A = get_molecule_info("Compound A", "N#CC(Nc1ccccc1)c1ncccc1O")
    
    # Print the identification of Compound A
    print("--- Identification of Compound A ---")
    print(f"The reaction produces an α-aminonitrile.")
    print(f"Product Name: {compound_A['name']} (2-(cyano(phenylamino)methyl)pyridin-3-ol)")
    print(f"SMILES String: {compound_A['smiles']}")
    print(f"Molecular Formula: {compound_A['formula']}")
    print(f"Molecular Weight: {compound_A['mw']:.2f} g/mol")
    print("\n" + "="*40 + "\n")

    # Print the overall chemical reaction
    # The overall reaction is: Reactant1 + Reactant2 + NaCN -> CompoundA + NaOH
    # To balance it, we consider the formation of H2O in the first step and the consumption
    # of a proton in the second, often simplified to: R1+R2+HCN -> Product
    # The problem description suggests TFE as a solvent, which can act as a proton source.
    # Water is formed in the first step (condensation).
    print("--- Overall Reaction Summary ---")
    
    # Balancing the atoms requires considering the full process.
    # Step 1: C6H5NO2 + C6H7N -> C12H10N2O (Imine) + H2O
    # Step 2: C12H10N2O + NaCN + H+ -> C13H11N3O (Compound A) + Na+
    # Overall: C6H5NO2 + C6H7N + NaCN + [H+] -> C13H11N3O + H2O + [Na+]
    # where [H+] is supplied by the acidic environment (TFE).
    
    print("Balanced Chemical Equation:")
    print(f"{reactant1['name']} + {reactant2['name']} + {reactant3_name} --(TFE)--> {compound_A['name']} + H2O + Na+")
    print(f"{reactant1['formula']} + {reactant2['formula']} + NaCN --> {compound_A['formula']} + H2O + Na+")
    
    # As requested, outputting the numbers (molecular weights) in the equation
    print("\nEquation with Molecular Weights (g/mol):")
    reactant3_mw = 22.99 + 12.01 + 14.01 # NaCN
    water_mw = 18.02
    na_ion_mw = 22.99
    
    # The final equation line showing all numbers
    print(f"{reactant1['mw']:.2f} + {reactant2['mw']:.2f} + {reactant3_mw:.2f}  -->  {compound_A['mw']:.2f} + {water_mw:.2f} + {na_ion_mw:.2f} (excluding proton source)")
    

if __name__ == "__main__":
    main()
