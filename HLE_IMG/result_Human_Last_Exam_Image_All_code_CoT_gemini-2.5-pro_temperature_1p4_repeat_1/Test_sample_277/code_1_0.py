def classify_molecules():
    """
    An interactive guide to classify the relationship between two molecules.
    Please answer the questions with 'yes' or 'no'.
    """
    print("--- Molecule Relationship Classifier ---")
    print("Please examine the two molecules and answer the following questions to find the correct choice.")
    print("The options are: (a) conformers isomers, (b) constitutional isomers, (c) Identical, (d) stereoisomers, (e) None of these\n")

    # Step 1: Compare Molecular Formula
    formula_same = input("1. Do the two molecules have the same molecular formula (same count of each type of atom)? (yes/no): ").strip().lower()
    if formula_same != 'yes':
        print("\nResult: The molecules are not isomers because their molecular formulas are different.")
        print("<<<e>>>")
        return

    # Step 2: Compare Connectivity
    connectivity_same = input("2. Do the two molecules have the same connectivity (are the atoms bonded together in the same order)? (yes/no): ").strip().lower()
    if connectivity_same != 'yes':
        print("\nResult: The molecules are constitutional isomers. They have the same formula but different bond connectivity.")
        print("<<<b>>>")
        return

    # At this point, formula and connectivity are the same.
    print("\nSince the formula and connectivity are the same, we must distinguish between identical molecules and stereoisomers.")

    # Step 3: Check for superimposability
    superimposable = input("3. Are the two molecules superimposable (can you place one perfectly on top of the other, perhaps after rotating it)? (yes/no): ").strip().lower()
    if superimposable == 'yes':
        print("\nResult: The molecules are identical. They are just different views of the same molecule.")
        print("<<<c>>>")
        return

    # At this point, they are non-superimposable stereoisomers.
    print("\nSince the molecules are not superimposable, they are a type of stereoisomer.")

    # Step 4: Check for interconversion by bond rotation
    bond_rotation = input("4. Can you convert one molecule into the other just by rotating around single bonds (without breaking any bonds)? (yes/no): ").strip().lower()
    if bond_rotation == 'yes':
        print("\nResult: The molecules are conformers (conformational isomers). They are different rotational states of the same molecule.")
        print("<<<a>>>")
    else:
        print("\nResult: The molecules are stereoisomers (specifically, configurational isomers). They have a different 3D arrangement that cannot be changed by simple bond rotation.")
        print("<<<d>>>")

if __name__ == '__main__':
    classify_molecules()