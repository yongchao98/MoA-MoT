# Import the necessary libraries from the RDKit package for chemical analysis
from rdkit import Chem
from rdkit.Chem import Descriptors
import re

def analyze_reaction():
    """
    Analyzes the starting material and predicts the product of the given reaction.
    """
    # The starting material is defined by its SMILES string.
    smiles_start = 'C[C](C)(C)[Si](C)(C)O[C@@H]1C/C=C(\\C[C@@H]1)/[C@]2(O)[C@H]3[C@@H](C=C[C@H]2C3)C(OC)OC'
    reactant_name = "(1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol"
    
    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_start)

    print("--- Reaction Analysis ---")
    print(f"Starting Material: {reactant_name}")

    # --- Step 1: Analyze functional groups and reactive motifs ---
    # Define a SMARTS pattern for a 3-hydroxy-1,5-diene, the substrate for an oxy-Cope rearrangement.
    # Pattern: [C=C]-[C with an OH group]-[any C]-[C=C]
    oxy_cope_pattern = Chem.MolFromSmarts('[#6]=[#6]-[#6X4]([#8H1])-[#6X4]-[#6]=[#6]')
    
    if mol.HasSubstructMatch(oxy_cope_pattern):
        print("\nObservation: The starting material contains a 3-hydroxy-1,5-diene system.")
        print("This indicates the reaction is an anionic oxy-Cope rearrangement.")
    else:
        print("\nWarning: The specific substructure for an oxy-Cope rearrangement was not detected automatically.")
        # Fallback explanation if pattern fails
        print("However, based on the reagents (KH), the tertiary alcohol adjacent to two alkene systems will undergo an anionic oxy-Cope rearrangement.")

    # --- Step 2: Describe the reaction mechanism ---
    print("\nReaction Mechanism Steps:")
    print("1. Deprotonation: KH, a strong base, deprotonates the tertiary alcohol (-OH) to form a potassium alkoxide (-O-).")
    print("2. [3,3]-Sigmatropic Rearrangement: The alkoxide undergoes a rapid, concerted rearrangement. This ring-expanding isomerization cleaves a C-C bond in the norbornene skeleton and forms a new C-C bond, creating a larger fused-ring system and an enolate intermediate.")
    print("3. Tautomerization: The aqueous workup (H2O/MeOH) protonates the enolate, which tautomerizes to the stable ketone product.")

    # --- Step 3: Describe the product and unchanged groups ---
    print("\nProduct Description:")
    print("- The tertiary alcohol group is converted into a ketone (C=O).")
    print("- The carbon skeleton is expanded into a new, complex bicyclic system.")
    print("- The tert-butyldimethylsilyl (TBDMS) ether and the dimethoxy ketal are stable to these conditions and remain intact.")

    # --- Step 4: Present the final chemical equation with elemental count ---
    print("\n--- Final Equation ---")
    # This reaction is an isomerization, so the molecular formula does not change.
    mol_formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    product_description = "An isomeric bicyclic ketone formed via ring expansion"
    
    print(f"Reactant: {reactant_name}")
    print("     |")
    print("     V  [ 1. KH, THF; 2. H2O/MeOH ]")
    print(f"Product:  {product_description}")
    
    print("\nThe final equation shows the conservation of atoms (isomerization):")
    
    # Extract element and count for detailed printing
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', mol_formula)
    equation_parts = []
    for element, count in elements:
        # If no number follows an element, its count is 1
        if count == '':
            count = '1'
        equation_parts.append(f"{element} = {count}")
    
    print(f"Reactant Formula: {mol_formula}")
    print(f"Product Formula:  {mol_formula}")
    print("Elemental Composition of Reactant and Product:")
    # Print each number in the final "equation" as requested
    print(', '.join(equation_parts))

# Execute the analysis function
analyze_reaction()