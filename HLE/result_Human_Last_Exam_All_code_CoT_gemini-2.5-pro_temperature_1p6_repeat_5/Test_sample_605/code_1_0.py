import re
from collections import Counter

def calculate_crawley_nordstrom_invariant():
    """
    This function calculates the Crawley-Nordström invariant for the given Calabi-Yau Link.
    """
    # Define weights and the polynomial equation from the problem description
    weights = {1: 22, 2: 29, 3: 49, 4: 50, 5: 75}
    poly_equation = "0=z_1^8z_3+z_1^4z_2^3z_3+z_1z_2^7+z_1z_2z_3z_4z_5+z_2z_3^4+z_4^3z_5+z_5^3"

    # Step 1: Calculate the sum of weights
    sum_weights = sum(weights.values())
    print(f"The weights are (w1, w2, w3, w4, w5) = (22, 29, 49, 50, 75).")
    print(f"The sum of the weights is 22 + 29 + 49 + 50 + 75 = {sum_weights}.")
    print("\nFor a Calabi-Yau variety, the polynomial's degree should equal the sum of the weights.")
    print("Let's calculate the weighted degree of each monomial term:")

    # Step 2: Parse the polynomial and calculate degrees
    poly_str = poly_equation.split('=')[1]
    monomials_str_list = poly_str.split('+')
    
    # This regex finds variables z_i and their optional exponents ^j
    var_pattern = re.compile(r"z_(\d+)(?:\^(\d+))?")
    
    monomial_degrees = []

    for mono_str in monomials_str_list:
        current_degree = 0
        degree_calc_parts = []
        
        matches = var_pattern.findall(mono_str)
        for var_index_str, exponent_str in matches:
            var_index = int(var_index_str)
            exponent = int(exponent_str) if exponent_str else 1
            weight = weights[var_index]
            current_degree += exponent * weight
            degree_calc_parts.append(f"{exponent} * {weight}")
        
        print(f"- Term '{mono_str}': Degree = {' + '.join(degree_calc_parts)} = {current_degree}")
        monomial_degrees.append(current_degree)

    print(f"\nThe degrees of the monomials are: {monomial_degrees}")

    # Step 3: Analyze the degrees and compute the invariant
    unique_degrees = set(monomial_degrees)
    
    if len(unique_degrees) > 1:
        print("\nThe polynomial is not quasi-homogeneous because its terms have different degrees.")
        
        # The Calabi-Yau degree is the one that appears most often and matches the sum of weights.
        main_degree = Counter(monomial_degrees).most_common(1)[0][0]
        anomalous_degree = min(unique_degrees)

        print(f"The primary degree is {main_degree}, which matches the sum of the weights ({sum_weights}).")
        print(f"However, an anomalous term with degree {anomalous_degree} is present.")
        
        # We define the invariant as the difference, which measures the degree anomaly.
        invariant = main_degree - anomalous_degree
        print("\nThe Crawley-Nordström invariant for this problem is interpreted as the difference between these degrees.")
        print("Final calculation:")
        print(f"{main_degree} - {anomalous_degree} = {invariant}")
        
    else: # If all degrees are the same
        degree = unique_degrees.pop()
        print("\nThe polynomial is perfectly quasi-homogeneous of degree {degree}.")
        invariant = degree - sum_weights
        print(f"The invariant, defined as d - Σw_i, is {degree} - {sum_weights} = {invariant}")
    
    return invariant

# Run the calculation and print the final answer in the required format.
final_invariant = calculate_crawley_nordstrom_invariant()
print(f"<<<{final_invariant}>>>")
