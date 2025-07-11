def print_obstruction_groups():
    """
    This function prints the list of homology and homotopy groups
    that serve as the building blocks for the homotopy-theoretic
    obstructions mentioned in the problem.
    """
    n = "n"
    k = "k"
    q = "q"

    print("The homotopy-theoretic obstructions are elements of the group pi_1(Aut(E)).")
    print("This group's structure, and thus the classification of obstructions, is determined by the following list of groups:")
    print("-" * 70)
    
    # Obstruction Group 1
    group1 = f"pi_1(SO(2*{k}))"
    description1 = "This is the first homotopy group of the special orthogonal group. The primary obstruction lies in a subgroup of this group."
    print(f"1. {group1}")
    print(f"   Description: {description1}")
    print("")

    # Obstruction Group 2
    group2 = f"pi_2(SO(2*{k}))"
    description2 = "This is the second homotopy group of the special orthogonal group. It maps into the space of higher obstructions."
    print(f"2. {group2}")
    print(f"   Description: {description2}")
    print("")

    # Obstruction Group 3
    group3 = f"H_tilde^{int(n)-1}(X; pi_{q}(SO(2*{k}))) for {q} >= 1"
    description3 = f"These are the reduced cohomology groups of the space X. Since X is a homology ({n}-1)-sphere, this is the only non-trivial cohomology dimension. The coefficients are the homotopy groups of SO(2*{k}). These groups determine the higher-order obstructions."
    print(f"3. {group3}")
    print(f"   Description: {description3}")
    print("-" * 70)

print_obstruction_groups()