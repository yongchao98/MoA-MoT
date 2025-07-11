import collections

def solve_cobordism_elements():
    """
    Calculates the number of distinct homology cobordism elements
    representable by integral surgery on knots with at most four crossings.
    """
    
    # Step 1: Define knots with at most 4 crossings and their invariants.
    # Knots: 0_1 (unknot), 3_1 (trefoil), 4_1 (figure-eight).
    # The trefoil (3_1) is chiral, so we consider both right-handed and left-handed versions.
    # The unknot and figure-eight are amphichiral.
    # We need the Arf invariant and the Tau (τ) invariant for each.
    # Arf(-K) = Arf(K), τ(-K) = -τ(K).
    Knot = collections.namedtuple('Knot', ['name', 'Arf', 'tau'])
    knots = [
        Knot("0_1 (Unknot)", 0, 0),
        Knot("3_1 (Right-handed Trefoil)", 1, 1),
        Knot("-3_1 (Left-handed Trefoil)", 1, -1),
        Knot("4_1 (Figure-Eight Knot)", 0, 0),
    ]

    print("Analyzing homology cobordism elements from integral surgery on knots with at most four crossings.")
    print("-" * 80)
    print("The knots considered are the unknot (0_1), trefoil (3_1), and figure-eight (4_1).")
    print("Integral surgery on a knot K gives a homology sphere for surgery coefficients +1 and -1.")
    print("We compute the invariants (d, μ) for each resulting homology sphere.")
    print("  - Rokhlin invariant: μ = Arf(K) mod 2")
    print("  - d-invariant for +1 surgery: d = -2 * τ(K)")
    print("  - d-invariant for -1 surgery: d = 2 * τ(K)")
    print("-" * 80)
    
    # Use a set to store the unique invariant pairs (d, μ), which identify the elements.
    represented_elements = set()

    # Step 2-4: Iterate through knots, calculate invariants for +/- 1 surgery.
    for knot in knots:
        print(f"Knot: {knot.name}")
        print(f"  Properties: Arf = {knot.Arf}, τ = {knot.tau}")
        
        # +1 Surgery
        mu_plus_1 = knot.Arf % 2
        d_plus_1 = -2 * knot.tau
        element_plus_1 = (d_plus_1, mu_plus_1)
        represented_elements.add(element_plus_1)
        print(f"  +1 Surgery: μ = {knot.Arf} mod 2 = {mu_plus_1}")
        print(f"              d = -2 * {knot.tau} = {d_plus_1}")
        print(f"              Generated element (d, μ): {element_plus_1}")
        
        # -1 Surgery
        mu_minus_1 = knot.Arf % 2
        d_minus_1 = 2 * knot.tau
        element_minus_1 = (d_minus_1, mu_minus_1)
        represented_elements.add(element_minus_1)
        print(f"  -1 Surgery: μ = {knot.Arf} mod 2 = {mu_minus_1}")
        print(f"              d = 2 * {knot.tau} = {d_minus_1}")
        print(f"              Generated element (d, μ): {element_minus_1}")
        print()

    # Step 5: Count and display the unique elements.
    print("-" * 80)
    print("The unique (d, μ) invariant pairs found are:")
    # Sort for deterministic output
    for elem in sorted(list(represented_elements)):
        description = ""
        if elem == (0, 0):
            description = "(This is the trivial element, the class of S^3)"
        elif elem == (2, 1) or elem == (-2, 1):
            description = "(These are non-trivial elements of infinite order)"
        print(f"  Element: (d={elem[0]}, μ={elem[1]}) {description}")

    count = len(represented_elements)
    print("\n" + "=" * 80)
    print(f"The total number of distinct elements of the homology cobordism group that can be represented is {count}.")
    print("=" * 80)

solve_cobordism_elements()