def solve_braid_index():
    """
    Calculates the braid index of the three-twist knot (5_2) to determine
    the tightest upper bound from Vogel's algorithm.
    """
    knot_name = "Three-twist knot (5_2)"
    
    print(f"Step 1: Find a lower bound for the braid index of the {knot_name}.")
    print("We use the Jones polynomial and the Morton-Franks-Williams inequality.")
    print("-" * 30)

    # The Jones polynomial for 5_2 is V(t) = t^2 - t + 1 - t^-1 + t^-2
    jones_poly_str = "V(t) = t^2 - t^1 + t^0 - t^-1 + t^-2"
    min_deg = -2
    max_deg = 2
    
    print(f"The Jones polynomial for the {knot_name} is: {jones_poly_str}")
    print(f"The minimum degree of the polynomial is: {min_deg}")
    print(f"The maximum degree of the polynomial is: {max_deg}")
    
    span = max_deg - min_deg
    print(f"The span of the polynomial is calculated as: max_degree - min_degree")
    print(f"span = {max_deg} - ({min_deg}) = {span}")
    
    print("\nThe Morton-Franks-Williams inequality states: span(V) <= 2 * b(K) - 2")
    print(f"Substituting the span for the {knot_name}:")
    
    equation_val_1 = span
    print(f"{equation_val_1} <= 2 * b(5_2) - 2")
    
    equation_val_2 = span + 2
    print(f"{equation_val_1} + 2 <= 2 * b(5_2)")
    print(f"{equation_val_2} <= 2 * b(5_2)")

    lower_bound = equation_val_2 / 2
    print(f"{equation_val_2} / 2 <= b(5_2)")
    print(f"Therefore, the braid index b(5_2) must be greater than or equal to {int(lower_bound)}.")
    print("-" * 30)

    print("\nStep 2: Find an upper bound for the braid index.")
    print("This can be done by finding an explicit braid representation for the knot.")
    print("-" * 30)

    braid_representation = "σ₁² σ₂⁻¹ σ₁⁻¹ σ₂"
    num_strands = 3
    
    print(f"The {knot_name} has a known braid representation on {num_strands} strands.")
    print(f"One such representation is given by the braid word: {braid_representation}")
    print(f"The existence of a {num_strands}-strand representation proves that the braid index is at most {num_strands}.")
    print(f"So, b(5_2) <= {num_strands}.")
    print("-" * 30)

    print("\nStep 3: Combine the bounds to find the exact braid index.")
    print("-" * 30)
    print(f"From Step 1, we have the lower bound: b(5_2) >= {int(lower_bound)}")
    print(f"From Step 2, we have the upper bound: b(5_2) <= {num_strands}")
    
    final_braid_index = 3
    print(f"\nCombining these two results, the only integer value possible is {final_braid_index}.")
    print(f"So, the braid index of the {knot_name} is exactly {final_braid_index}.")
    
    print("\nConclusion:")
    print("Vogel's algorithm provides an upper bound for the braid index. The tightest possible")
    print(f"upper bound is the braid index itself, which we have found to be {final_braid_index}.")

solve_braid_index()