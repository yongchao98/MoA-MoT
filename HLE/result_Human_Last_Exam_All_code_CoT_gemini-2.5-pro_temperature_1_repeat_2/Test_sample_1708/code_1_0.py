def solve_order_type():
    """
    This function explains the derivation of the order type for the set of
    finite strings over {a, b, c, d} under lexicographical order.
    """

    # Let L be the order type we are looking for.
    # The set W is all finite strings from {a,b,c,d}, ordered lexicographically.

    # Step 1: Decompose the set W into four partitions.
    # S_a = strings starting with 'a'
    # S_b = strings starting with 'b'
    # S_c = strings starting with 'c'
    # S_d = strings starting with 'd'
    # The total order type L is the ordinal sum of the order types of these partitions.
    # L = OT(S_a) + OT(S_b) + OT(S_c) + OT(S_d)

    # Step 2: Find the order type of each partition.
    # Consider S_a. It contains "a", "aa", "ab", etc.
    # The structure of S_a is the single element "a" followed by an ordered set
    # that is isomorphic to the original set W (e.g., via the map s -> "as").
    # The order type of a set with one minimum element followed by a set of type L
    # is given by the ordinal expression 1 + L.
    # By symmetry, all four partitions have this same order type.
    # OT(S_a) = OT(S_b) = OT(S_c) = OT(S_d) = 1 + L

    # Step 3: Formulate the final equation for L.
    # Substituting this back into our sum gives us an equation for L.
    print("The order type, L, is defined by the following ordinal equation:")
    one = 1
    four = 4
    print(f"L = (1 + L) + (1 + L) + (1 + L) + (1 + L)")
    print(f"L = ({one} + L) * {four}")

    # Step 4: Solve the ordinal equation.
    # The set W is infinite and has no largest element (for any string s, "sd" is larger).
    # This means L must be a limit ordinal.
    # For any limit ordinal L, a key property is that 1 + L = L.
    # So, the equation simplifies.
    print("\nBecause L must be a limit ordinal, 1 + L simplifies to L.")
    print("The equation becomes:")
    print(f"L = L * {four}")

    # Step 5: Identify the ordinal L.
    # We need to find the smallest non-zero limit ordinal L that satisfies L = L * 4.
    # Let's test the first few limit ordinals:
    # - If L = ω (omega), then ω * 4 = ω+ω+ω+ω, which is strictly greater than ω.
    # - If L = ω^2, then ω^2 * 4 is strictly greater than ω^2.
    # The smallest non-zero ordinal α such that α = α * k for a finite k > 1
    # is the ordinal ω^ω (omega to the power of omega).

    solution = "ω^ω"
    print(f"\nThe smallest non-zero ordinal that solves the equation L = L * 4 is {solution}.")
    
solve_order_type()