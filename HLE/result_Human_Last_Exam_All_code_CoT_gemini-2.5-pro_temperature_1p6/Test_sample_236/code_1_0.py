def solve_homology_cobordism_elements():
    """
    Calculates how many homology cobordism group elements can be represented
    by an integral surgery on a knot with at most four crossings.
    """
    # A dictionary mapping knots (with at most 4 crossings) to their Arf invariants.
    # The Arf invariant determines the element in the homology cobordism group (Z/2Z).
    # Arf(K) = 0 -> Trivial element
    # Arf(K) = 1 -> Non-trivial element
    knots_and_arfs = {
        "0_1 (Unknot)": 0,
        "3_1 (Trefoil)": 1,
        "4_1 (Figure-eight)": 0
    }

    print("Analysis of Homology Cobordism Group Elements from Knots <= 4 Crossings\n")
    print("The element is determined by the knot's Arf invariant.")
    print("-" * 60)

    # Use a set to automatically handle uniqueness of the resulting group elements (0 or 1)
    found_elements = set()
    arf_values_for_equation = []

    for knot, arf in knots_and_arfs.items():
        found_elements.add(arf)
        arf_values_for_equation.append(str(arf))
        element_type = "Trivial" if arf == 0 else "Non-trivial"
        
        print(f"Knot: {knot}")
        print(f"Arf Invariant: {arf}")
        print(f"Resulting Element: {element_type} (class {arf})")
        print("-" * 60)
        
    num_elements = len(found_elements)
    
    # Construct the final equation string showing all numbers involved
    arfs_str = ", ".join(arf_values_for_equation)
    final_equation = f"count(unique values of {{{arfs_str}}}) = {num_elements}"

    print("Final Result:")
    print(f"The distinct elements found correspond to Arf invariants: {sorted(list(found_elements))}")
    print("The final count is derived from the unique Arf invariants:")
    print(final_equation)


solve_homology_cobordism_elements()
<<<2>>>