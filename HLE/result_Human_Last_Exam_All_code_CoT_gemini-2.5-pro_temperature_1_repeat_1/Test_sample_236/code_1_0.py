import sys

def solve():
    """
    Calculates how many elements of the homology cobordism group can be
    represented by integral surgery on knots with at most four crossings.
    """

    # Knots with at most 4 crossings and their standard Alexander polynomials.
    # The polynomial is represented as a lambda function of a variable 't'.
    knots = {
        "Unknot (0_1)": lambda t: 1,
        "Trefoil (3_1)": lambda t: t**2 - t + 1,
        "Figure-eight (4_1)": lambda t: t**2 - 3*t + 1,
    }

    print("Analyzing knots with at most four crossings to determine which elements of the")
    print("homology cobordism group (Theta^3_H) can be represented by integral surgery.")
    print("This group has two elements, distinguished by the Rokhlin invariant (0 or 1).\n")
    print("The Rokhlin invariant of a homology sphere from +/-1 surgery on a knot K is equal to the Arf invariant of K.")
    print("The Arf invariant is calculated from the knot's Alexander polynomial evaluated at t=-1.\n")

    represented_elements = set()

    for name, alexander_poly in knots.items():
        print(f"--- Analyzing Knot: {name} ---")

        # The Alexander polynomial is defined up to units +/- t^k.
        # We use the standard, symmetrized version.
        # For example, for the Trefoil, Delta(t) = t^2 - t + 1.
        # For the Figure-eight, Delta(t) = t^2 - 3t + 1.

        # Evaluate the Alexander polynomial at t = -1.
        # This value is related to the knot determinant.
        determinant = alexander_poly(-1)
        print(f"Alexander polynomial at t=-1 is: {determinant}")

        # Calculate the Arf invariant based on the determinant modulo 8.
        # Arf = 0 if det = +/- 1 (mod 8)
        # Arf = 1 if det = +/- 3 (mod 8)
        arf_invariant = -1 # Default/error value
        det_mod_8 = determinant % 8
        if det_mod_8 in (1, 7): # 7 is congruent to -1 mod 8
            arf_invariant = 0
        elif det_mod_8 in (3, 5): # 5 is congruent to -3 mod 8
            arf_invariant = 1

        if arf_invariant == 0:
            print(f"The determinant ({determinant}) is congruent to {determinant % 8} mod 8.")
            print(f"This implies the Arf invariant is 0.")
            print(f"Therefore, surgery on this knot represents the TRIVIAL element of Theta^3_H.")
            represented_elements.add("Trivial")
        elif arf_invariant == 1:
            print(f"The determinant ({determinant}) is congruent to {determinant % 8} mod 8.")
            print(f"This implies the Arf invariant is 1.")
            print(f"Therefore, surgery on this knot represents the NON-TRIVIAL element of Theta^3_H.")
            represented_elements.add("Non-Trivial")
        else:
            print("Could not determine the Arf invariant from the determinant.")
        
        print("-" * (len(name) + 22) + "\n")

    num_elements = len(represented_elements)

    print("Summary:")
    print(f"The set of representable elements is: {represented_elements}")
    print(f"The total number of distinct elements of the homology cobordism group that can be")
    print(f"represented by integral surgery on a knot with at most four crossings is: {num_elements}")

    # The final answer must be wrapped in <<<>>>
    # We use sys.stdout.write to avoid adding a newline character.
    sys.stdout.write(f"<<<{num_elements}>>>")

solve()