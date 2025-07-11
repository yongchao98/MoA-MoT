def check_equivalence():
    """
    This function demonstrates that the two minimal grid diagrams for the trefoil knot
    are equivalent under translation.
    """
    # The two permutations representing the minimal (3x3) trefoil knot.
    p1 = (2, 3, 1)
    p2 = (3, 1, 2)

    print(f"The set of minimal grid diagrams for the left-hand trefoil knot corresponds to two permutations:")
    print(f"Diagram 1: p1 = {p1}")
    print(f"Diagram 2: p2 = {p2}")
    print("\nWe check if they are equivalent under translation.")
    print("A translation can be a cyclic shift of the columns.")
    
    # Apply a vertical translation by d=1 to p1.
    # The formula is p'_i = (p_i - 1 + d) % n + 1
    # Here, n=3 and d=1.
    
    p1_translated = []
    d = 1
    n = 3
    for x in p1:
        # We need to perform modulo arithmetic on a 0-indexed base
        # then convert back to 1-indexed.
        translated_val = (x - 1 + d) % n + 1
        p1_translated.append(translated_val)
    
    p1_translated = tuple(p1_translated)

    print(f"\nApplying a column shift of d={d} to p1 = {p1}...")
    print(f"The calculation is new_p_i = (old_p_i - 1 + {d}) % {n} + 1 for each element.")
    
    # Show the calculation step-by-step
    new_p_1 = (p1[0] - 1 + d) % n + 1
    new_p_2 = (p1[1] - 1 + d) % n + 1
    new_p_3 = (p1[2] - 1 + d) % n + 1

    print(f"For the first element 2: ({p1[0]} - 1 + {d}) % {n} + 1 = {new_p_1}")
    print(f"For the second element 3: ({p1[1]} - 1 + {d}) % {n} + 1 = {new_p_2}")
    print(f"For the third element 1: ({p1[2]} - 1 + {d}) % {n} + 1 = {new_p_3}")

    print(f"The resulting permutation is {p1_translated}.")

    if p1_translated == p2:
        print(f"\nThis result, {p1_translated}, is identical to p2, {p2}.")
        print("This means the two diagrams are in the same equivalence class.")
        print("Since these are the only two minimal diagrams, there is only 1 unique diagram up to translation.")
    else:
        print("\nThe diagrams are not equivalent by this translation.")

check_equivalence()