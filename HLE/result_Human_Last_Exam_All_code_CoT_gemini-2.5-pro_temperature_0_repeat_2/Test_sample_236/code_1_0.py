def count_homology_cobordism_elements():
    """
    Calculates the number of homology cobordism group elements representable
    by an integral surgery on a knot with at most four crossings.
    """

    # Step 1: Define knots with at most 4 crossings and their signatures.
    # We consider right-handed (RH) and left-handed (LH) versions for chiral knots.
    # The figure-eight knot (4_1) is amphichiral, so one entry suffices.
    knots = {
        'Unknot (0_1)': 0,
        'Right-handed Trefoil (3_1)': -2,
        'Left-handed Trefoil (-3_1)': 2,
        'Figure-eight knot (4_1)': 0
    }

    # This set will store the unique elements found. We use strings for representation.
    unique_elements = set()

    print("Analysis of Homology Cobordism Elements from Knot Surgery")
    print("=" * 60)
    print("A surgery S^3_n(K) is a homology sphere if n=±1.")
    print("The resulting element is trivial if the knot signature σ(K) is a multiple of 8.")
    print("-" * 60)

    # Step 2 & 3: Iterate through knots and determine the resulting elements.
    for name, signature in knots.items():
        # Check if the signature is a multiple of 8.
        if signature % 8 == 0:
            # If so, both +1 and -1 surgeries give the trivial element.
            element = "Identity Element"
            unique_elements.add(element)
            print(f"Knot: {name:28} | σ = {signature:2} | Result: {element}")
        else:
            # Otherwise, we get non-trivial elements.
            # For the trefoil (σ=±2), we get the Poincaré sphere (P) and its inverse (-P).
            # We use a canonical representation for these.
            element1 = "Poincaré Sphere (P)"
            element2 = "Inverse Poincaré Sphere (-P)"
            unique_elements.add(element1)
            unique_elements.add(element2)
            print(f"Knot: {name:28} | σ = {signature:2} | Results: {element1} and {element2}")

    # Step 4: Count the unique elements.
    count = len(unique_elements)

    print("=" * 60)
    print("The distinct elements found are:")
    for element in sorted(list(unique_elements)):
        print(f"- {element}")

    print("\nFinal Count:")
    # The problem asks to output the numbers in the final equation.
    # Here, the "equation" is the summation of the unique elements found.
    # We have 1 (Identity) + 1 (P) + 1 (-P) = 3
    print("1 (from Unknot/Fig-8) + 2 (from Trefoil) = 3 distinct elements.")
    print(f"Total number of elements: {count}")

if __name__ == '__main__':
    count_homology_cobordism_elements()