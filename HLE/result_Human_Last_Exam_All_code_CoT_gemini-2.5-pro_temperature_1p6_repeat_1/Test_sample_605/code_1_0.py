import collections

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a given Calabi-Yau link.
    """
    # Step 1: Define the weights and the exponents of the monomials from the polynomial
    # 0 = z_1^8*z_3 + z_1^4*z_2^3*z_3 + z_1*z_2^7 + z_1*z_2*z_3*z_4*z_5 + z_2*z_3^4 + z_4^3*z_5 + z_5^3
    weights = [22, 29, 49, 50, 75]
    exponents_list = [
        [8, 0, 1, 0, 0],  # z_1^8*z_3
        [4, 3, 1, 0, 0],  # z_1^4*z_2^3*z_3
        [1, 7, 0, 0, 0],  # z_1*z_2^7
        [1, 1, 1, 1, 1],  # z_1*z_2*z_3*z_4*z_5
        [0, 1, 4, 0, 0],  # z_2*z_3^4
        [0, 0, 0, 3, 1],  # z_4^3*z_5
        [0, 0, 0, 0, 3]   # z_5^3
    ]

    # Step 2: Calculate the weighted degree for each monomial
    degrees = []
    for exponents in exponents_list:
        degree = sum(e * w for e, w in zip(exponents, weights))
        degrees.append(degree)

    # For the polynomial to be well-defined, it must be weighted-homogeneous.
    # We assume the degree 'd' is the most common degree among the terms,
    # accounting for potential typos in the input polynomial.
    if not degrees:
        d = 0
    else:
        degree_counts = collections.Counter(degrees)
        d = degree_counts.most_common(1)[0][0]

    # Step 3: Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # Step 4: Calculate the Crawley-Nordstrom invariant
    invariant = d - sum_of_weights

    # Output the result in the specified equation format
    weight_sum_str = " + ".join(map(str, weights))
    print(f"The degree of the polynomial is d = {d}.")
    print(f"The sum of the weights is S = {weight_sum_str} = {sum_of_weights}.")
    print("The Crawley-Nordström invariant is calculated as c = d - S.")
    print("The final equation is:")
    print(f"{d} - ({weight_sum_str}) = {invariant}")

    # Return the final numerical answer
    return invariant

# Run the calculation and store the result
result = calculate_crawley_nordstrom_invariant()
# Print the final answer in the required format
print(f"<<<{result}>>>")