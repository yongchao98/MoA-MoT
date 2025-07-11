import collections

def solve_molecule_rotation():
    """
    Analyzes the molecular structure to determine which phenyl rings can rotate freely.
    """
    # Step 1: Represent the phenyl rings and their structural role.
    # 'fixed' means it is part of a rigid cyclometalated chelate ring.
    # 'peripheral' means it is a substituent attached by a single bond.
    phenyl_rings = {
        'A': 'fixed',
        'C': 'fixed',
        'G': 'peripheral',
        'H': 'peripheral',
        'J': 'peripheral',
        'K': 'peripheral',
        'L': 'peripheral'
    }

    print("Analyzing which phenyl rings can rotate freely in solution.")
    print("Principle: Rotation is free around single C-C bonds unless the ring is part of a rigid structure (e.g., a chelate ring).\n")

    freely_rotating_rings = []
    
    # Using an ordered dictionary to process in a defined order for clear output.
    # This is not strictly necessary for the logic but makes the printout cleaner.
    ordered_rings = collections.OrderedDict(sorted(phenyl_rings.items()))

    # Step 2: Iterate through each phenyl ring and apply the principle.
    for ring, role in ordered_rings.items():
        if role == 'fixed':
            print(f"Ring {ring}: This ring is cyclometalated and part of a rigid chelate system. It CANNOT rotate freely.")
        elif role == 'peripheral':
            print(f"Ring {ring}: This ring is a peripheral substituent attached by a single C-C bond. It CAN rotate freely.")
            freely_rotating_rings.append(ring)

    # Step 3: Output the final list of rings.
    print("\n-----------------------------------------------------")
    print("Final Result:")
    print("The set of phenyl rings that can rotate freely is composed of the following individual rings:")
    # "output each number in the final equation!" is interpreted as outputting each label.
    for ring in freely_rotating_rings:
        print(f"- Ring {ring}")
    
    # Matching the result to the provided options
    # The identified rings are G, H, J, K, L.
    # This corresponds to option A: G, J, L, K, H (same set of rings).

solve_molecule_rotation()