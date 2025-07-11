def analyze_complex_lifetimes():
    """
    Analyzes the chemical structures of four iridium complexes to predict which
    will have shorter operational lifetimes in light-emitting electrochemical cells (LECs).
    """

    # Define the key feature of each complex: the number of fluorine atoms on each cyclometalating ligand.
    complex_features = {
        1: 2,  # di-fluorinated phenylpyridine ligand
        2: 2,  # di-fluorinated phenylpyridine ligand
        3: 1,  # mono-fluorinated phenylpyridine ligand
        4: 1   # mono-fluorinated phenylpyridine ligand
    }

    print("--- Analysis of Iridium Complex Stability and Lifetime ---")
    print("\nThe operational lifetime of an emitter in an LEC is largely determined by its chemical stability.")
    print("\nA key factor for the stability of these iridium complexes is the fluorination of the phenylpyridine ligands.")
    print("\nPrinciple: Fluorine atoms are strongly electron-withdrawing. A higher degree of fluorination increases the stability of the complex by:")
    print("1. Making it more resistant to oxidation (lowering HOMO energy).")
    print("2. Strengthening the chemical bonds within the ligand.")
    print("Therefore, more fluorine atoms lead to a longer lifetime, and fewer fluorine atoms lead to a shorter lifetime.")

    print("\n--- Comparing the Complexes ---")
    for complex_id, num_fluorines in complex_features.items():
        print(f"Complex {complex_id}: Has {num_fluorines} fluorine atom(s) per cyclometalating ligand.")

    # Identify the complexes with the least fluorination
    shorter_lifetime_complexes = []
    min_fluorines = min(complex_features.values())
    for complex_id, num_fluorines in complex_features.items():
        if num_fluorines == min_fluorines:
            shorter_lifetime_complexes.append(complex_id)

    print("\n--- Conclusion ---")
    print("Complexes with the lowest number of fluorine atoms are expected to be the least stable and have shorter lifetimes.")
    print("The complexes with shorter lifetimes are therefore:")
    # Print each number as requested
    print(shorter_lifetime_complexes[0])
    print(shorter_lifetime_complexes[1])

analyze_complex_lifetimes()