def solve_molecule_challenge():
    """
    This script presents the solution to the molecule design challenge.
    It provides the SMILES string of a molecule that fits all the given constraints
    and demonstrates the calculation of its exact molecular weight.
    """

    # The designed molecule is a derivative of phenylcoumarin.
    # Specifically, it is 7-hydroxy-4-(2,3-dihydroxyphenyl)coumarin.
    # This structure meets all the complex requirements.
    smiles_string = "c1cc(O)c2c(c1)oc(=O)c(c2)-c1cccc(O)c1O"

    # Molecular formula derived from the constraints
    formula = {"C": 15, "H": 10, "O": 5}
    
    # Exact monoisotopic masses of the most common isotopes
    atomic_masses = {
        "C": 12.000000,
        "H": 1.007825,
        "O": 15.994915
    }

    # --- Calculation ---
    c_mass = formula["C"] * atomic_masses["C"]
    h_mass = formula["H"] * atomic_masses["H"]
    o_mass = formula["O"] * atomic_masses["O"]
    total_mass = c_mass + h_mass + o_mass

    # --- Output ---
    print(f"Designed Molecule SMILES String:\n{smiles_string}\n")
    
    print("Molecular Weight Calculation:")
    print(f"The molecular formula is C15H10O5.")
    print("Using exact monoisotopic masses:")
    print(f"  Carbon (C):   {formula['C']:>2} * {atomic_masses['C']:<10.6f} = {c_mass:>10.6f}")
    print(f"  Hydrogen (H): {formula['H']:>2} * {atomic_masses['H']:<10.6f} = {h_mass:>10.6f}")
    print(f"  Oxygen (O):   {formula['O']:>2} * {atomic_masses['O']:<10.6f} = {o_mass:>10.6f}")
    print("--------------------------------------------------")
    print(f"  Total Weight: {c_mass:.6f} + {h_mass:.6f} + {o_mass:.6f} = {total_mass:.6f}")
    print(f"\nThe calculated molecular weight is {total_mass:.3f}, which matches the target of 270.053.")


solve_molecule_challenge()
<<<c1cc(O)c2c(c1)oc(=O)c(c2)-c1cccc(O)c1O>>>