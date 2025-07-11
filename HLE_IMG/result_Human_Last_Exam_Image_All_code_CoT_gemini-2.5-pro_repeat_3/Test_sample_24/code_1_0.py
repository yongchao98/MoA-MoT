import sys
try:
    from sympy.combinatorics.fp_groups import FpGroup
    from sympy.abc import x, y
except ImportError:
    print("Sympy library not found. Please install it using: pip install sympy")
    sys.exit(1)

def identify_manifold():
    """
    Identifies the 3-manifold from its fundamental group presentation derived
    from its Heegaard diagram.
    """
    # The Heegaard diagram is a known representation of the Poincaré homology sphere.
    # Its fundamental group is the binary icosahedral group. A presentation for this
    # group is G = <x, y | (xy)^2 = x^3, y^5 = x^3>.
    # We define the group by its relators (relations rewritten to be equal to 1).
    # R1: (xy)^2 * x^-3 = 1  => xyxyx^-3 = 1
    # R2: y^5 * x^-3 = 1      => y^5x^-3 = 1
    relators = [(x*y)**2 * x**-3, y**5 * x**-3]

    try:
        G = FpGroup(relators)
    except Exception as e:
        print(f"An error occurred while creating the group: {e}")
        return

    # The first homology group H1(M) is the abelianization of the fundamental group G.
    # Its structure is given by the abelian invariants.
    try:
        h1_invariants = G.abelian_invariants()
    except Exception as e:
        print(f"Could not compute abelian invariants: {e}")
        h1_invariants = "Computation failed"

    # We also compute the order of the fundamental group.
    # This can be computationally expensive but is feasible for this group.
    try:
        group_order = G.order()
    except Exception as e:
        print(f"Could not compute group order. It might be infinite or too complex: {e}")
        group_order = "Computation failed or group is infinite"


    print("The Heegaard diagram represents the Poincaré homology sphere, also known as the Brieskorn sphere Σ(2,3,5).")
    print("Its fundamental group, π1(M), is the binary icosahedral group.")
    print("A known presentation for this group is G = <x, y | R1, R2>, where the relations are:")
    # The final instruction "output each number in the final equation" is interpreted
    # as printing the numbers in the relations defining the group.
    print("Equation 1: (x*y)^2 = x^3")
    print("Equation 2: y^5 = x^3")

    print("\n--- Computational Verification ---")

    print("\n1. First Homology Group H1(M):")
    print("H1(M) is the abelianization of π1(M). Its structure is given by the abelian invariants.")
    print(f"The computed abelian invariants are: {h1_invariants}")
    if h1_invariants == [0] or h1_invariants == []:
        print("The trivial invariants ([0] or []) mean H1(M) is the trivial group {0}.")
        print("This confirms the manifold is a homology sphere.")
    else:
        print("The computed homology group is not trivial.")

    print("\n2. Fundamental Group π1(M):")
    print("We compute the order of the group to check if it is trivial (like S^3) or not.")
    print(f"The computed order of the fundamental group is: {group_order}")
    if group_order == 1:
        print("The fundamental group is trivial, suggesting the manifold is the 3-sphere S^3.")
    elif isinstance(group_order, int) and group_order > 1:
        print("The fundamental group is non-trivial, so the manifold is not the 3-sphere S^3.")

    print("\n--- Numbers in the defining equations ---")
    print("Equation 1, (x*y)^2 = x^3, contains the exponents: 2, 3")
    print("Equation 2, y^5 = x^3, contains the exponents: 5, 3")


if __name__ == "__main__":
    identify_manifold()