import collections

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.
    """
    # Step 1: Define weights and polynomial terms
    weights = [22, 29, 49, 50, 75]
    # Polynomial: z_1^8z_3 + z_1^4z_2^3z_3 + z_1z_2^7 + z_1z_2z_3z_4z_5 + z_2z_3^4 + z_4^3z_5 + z_5^3
    # Each inner list represents the exponents [a1, a2, a3, a4, a5] for a term.
    terms_exponents = [
        [8, 0, 1, 0, 0],  # z_1^8*z_3
        [4, 3, 1, 0, 0],  # z_1^4*z_2^3*z_3
        [1, 7, 0, 0, 0],  # z_1*z_2^7
        [1, 1, 1, 1, 1],  # z_1*z_2*z_3*z_4*z_5
        [0, 1, 4, 0, 0],  # z_2*z_3^4
        [0, 0, 0, 3, 1],  # z_4^3*z_5
        [0, 0, 0, 0, 3],  # z_5^3
    ]
    term_strings = [
        "z_1^8*z_3",
        "z_1^4*z_2^3*z_3",
        "z_1*z_2^7",
        "z_1*z_2*z_3*z_4*z_5",
        "z_2*z_3^4",
        "z_4^3*z_5",
        "z_5^3"
    ]

    print("Step 1: Verify quasi-homogeneity and find the polynomial's degree (d).")
    print(f"Weights (w1, w2, w3, w4, w5) = {tuple(weights)}")
    
    degrees = []
    # Calculate the weighted degree for each term
    for i, exponents in enumerate(terms_exponents):
        degree = sum(e * w for e, w in zip(exponents, weights))
        degrees.append(degree)
        calculation_parts = [f"{exp}*{w}" for exp, w in zip(exponents, weights) if exp > 0]
        calculation_str = " + ".join(calculation_parts)
        print(f"Degree of term '{term_strings[i]}': {calculation_str} = {degree}")

    # Determine the polynomial's degree d
    # Check if all degrees are the same
    degree_counts = collections.Counter(degrees)
    d, count = degree_counts.most_common(1)[0]
    
    if len(degree_counts) > 1:
        print(f"\nWarning: The polynomial is not quasi-homogeneous as given.")
        print(f"Most terms ({count}/{len(degrees)}) have a degree of {d}. We will proceed assuming d = {d}.")

    print(f"\nThe determined quasi-homogeneous degree of the polynomial is d = {d}.")
    
    # Step 2: Calculate the sum of weights
    sum_w = sum(weights)
    sum_w_str = " + ".join(map(str, weights))
    print(f"\nStep 2: Calculate the sum of the weights (Σw_i).")
    print(f"Σw_i = {sum_w_str} = {sum_w}.")

    # Step 3: Check the Calabi-Yau condition
    print("\nStep 3: Check the Calabi-Yau condition (d = Σw_i).")
    if d == sum_w:
        print(f"The condition is satisfied, as d = {d} and Σw_i = {sum_w}.")
    else:
        print(f"The condition is NOT satisfied, as d = {d} and Σw_i = {sum_w}.")

    # Step 4: Calculate the Crawley-Nordström invariant
    print("\nStep 4: Calculate the Crawley-Nordström invariant (cn).")
    print("Formula: cn(f) = (d - Σw_i) / d")
    
    numerator = d - sum_w
    
    # Print the final equation with values
    print(f"cn(f) = ({d} - {sum_w}) / {d}")
    print(f"cn(f) = {numerator} / {d}")
    
    if d == 0:
        cn = float('inf') if numerator > 0 else (float('-inf') if numerator < 0 else 0)
    else:
        cn = numerator / d
    
    print(f"cn(f) = {cn}")
    
    print(f"\n<<<0.0>>>")

calculate_crawley_nordstrom_invariant()