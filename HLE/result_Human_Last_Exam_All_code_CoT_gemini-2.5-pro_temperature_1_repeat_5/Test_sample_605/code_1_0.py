import collections

def solve_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau link.
    """
    # Step 1: Define weights and polynomial terms
    weights = [22, 29, 49, 50, 75]
    # Exponents for each term z_1^a1 * z_2^a2 * ...
    # f = z_1^8*z_3 + z_1^4*z_2^3*z_3 + z_1*z_2^7 + z_1*z_2*z_3*z_4*z_5 + z_2*z_3^4 + z_4^3*z_5 + z_5^3
    terms_exponents = [
        [8, 0, 1, 0, 0],
        [4, 3, 1, 0, 0],
        [1, 7, 0, 0, 0],
        [1, 1, 1, 1, 1],
        [0, 1, 4, 0, 0],
        [0, 0, 0, 3, 1],
        [0, 0, 0, 0, 3]
    ]

    # Step 2: Determine the degree d of the polynomial
    term_degrees = [sum(e * w for e, w in zip(term, weights)) for term in terms_exponents]

    if len(set(term_degrees)) > 1:
        print(f"Warning: The polynomial is not quasi-homogeneous. Calculated term degrees: {term_degrees}")
        # Assume the intended degree is the most common one
        degree_counts = collections.Counter(term_degrees)
        d = degree_counts.most_common(1)[0][0]
        print(f"Proceeding with the most common degree: d = {d}\n")
    else:
        d = term_degrees[0]
        print(f"The polynomial is quasi-homogeneous with degree d = {d}\n")

    # Step 3: Calculate the sum of the weights
    sum_of_weights = sum(weights)
    print(f"The sum of the weights is: Σwᵢ = {sum_of_weights}")
    print(f"The degree of the polynomial is: d = {d}")

    # Check the Calabi-Yau condition
    if d == sum_of_weights:
        print("The Calabi-Yau condition (d = Σwᵢ) is satisfied.\n")
    else:
        print("Warning: The Calabi-Yau condition (d = Σwᵢ) is not met.\n")

    # Step 4: Calculate the Crawley-Nordström invariant
    if d == 0:
        print("Error: Degree d is zero, cannot calculate the invariant.")
        return

    invariant = (d - sum_of_weights) / d

    print("The Crawley-Nordström invariant is calculated using the formula: (d - Σwᵢ) / d")
    # Using integer casting for cleaner output as the result is an integer
    print(f"Final Equation: ({d} - {sum_of_weights}) / {d} = {int(invariant)}")

solve_crawley_nordstrom_invariant()
<<<0>>>