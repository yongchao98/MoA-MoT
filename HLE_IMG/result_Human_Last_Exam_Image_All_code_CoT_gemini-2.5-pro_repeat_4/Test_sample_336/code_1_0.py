# Import necessary libraries from RDKit
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import GraphDescriptors

def solve_molecular_properties():
    """
    This function identifies molecules with logP > 1 from a given list,
    calculates their Eccentric Connectivity Index (including hydrogens),
    and prints the sum.
    """
    # Step 1: Define the molecules using their SMILES strings
    # The SMILES strings correspond to the three molecules in the image.
    smiles_and_names = {
        "Reactant 1": "Nc1nccc(F)c1OCc1cc(OC)c(OC)c(OC)c1",
        "Reactant 2": "ClN1CC2=C(C1)N=CN=C2",
        "Product": "Fc1cncc(Nc2ncc3c(n2)CC[NH]3)c1"
    }

    print("Analyzing molecules based on Crippen logP > 1...")
    
    eci_values_to_sum = []
    
    for name, smiles in smiles_and_names.items():
        # Create a molecule object from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse SMILES for {name}: {smiles}")
            continue

        # Step 2: Calculate Crippen logP
        logp = Crippen.MolLogP(mol)
        
        print(f"\n- Molecule: {name}")
        print(f"  SMILES: {smiles}")
        print(f"  Crippen logP = {logp:.4f}")

        # Step 3: Filter molecules with logP > 1
        if logp > 1:
            print("  Condition logP > 1 is met. Calculating ECI.")
            
            # Step 4: Calculate Eccentric Connectivity Index (ECI)
            # Add explicit hydrogens to the molecular graph
            mol_with_hs = Chem.AddHs(mol)
            
            # Calculate ECI on the graph with hydrogens
            eci = GraphDescriptors.EccentricConnectivityIndex(mol_with_hs)
            eci_values_to_sum.append(eci)
            
            print(f"  Eccentric Connectivity Index (with H's) = {eci}")
        else:
            print("  Condition logP > 1 is NOT met. Excluding from sum.")

    # Step 5: Sum the ECI values and print the result
    print("\n--------------------------------------------------")
    print("Final Calculation:")
    
    if not eci_values_to_sum:
        print("No molecules met the logP > 1 criteria.")
        total_eci = 0
        print("Sum of Eccentric Connectivity Indices = 0")
    else:
        # Format the equation string showing each number
        equation_str = " + ".join(map(str, [int(v) for v in eci_values_to_sum]))
        total_eci = sum(eci_values_to_sum)
        
        print(f"Sum of Eccentric Connectivity Indices = {equation_str}")
        print(f"Total Sum = {int(total_eci)}")

# Execute the function
solve_molecular_properties()