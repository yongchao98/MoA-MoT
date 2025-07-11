import math

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau link.
    """
    # Step 1: Define weights and degree
    # The weights of the ambient space C^5
    weights = [22, 29, 49, 50, 75]
    n = len(weights)
    
    # The degree 'd' is the sum of the weights for a Calabi-Yau hypersurface.
    # This is consistent with most monomials in the user-provided polynomial.
    d = sum(weights)
    
    print("This script calculates the Crawley-Nordström invariant for a Calabi-Yau link.")
    print(f"The calculation assumes a weighted-homogeneous polynomial in C^{n} with:")
    print(f"Weights w = {weights}")
    print(f"Degree  d = sum(w_i) = {d}\n")
    
    # Step 2: Calculate the parameters b_i = w_i / gcd(d, w_i)
    g = [math.gcd(d, w) for w in weights]
    b = [w // gi for w, gi in zip(weights, g)]
    
    print("The formula for the invariant requires the following parameters:")
    print(f"g_i = gcd(d, w_i): {g}")
    print(f"b_i = w_i / g_i:   {b}\n")
    
    # Step 3: Calculate the invariant using the combinatorial formula.
    # The invariant is Sum_{j>=0} Product_{i=1..n} C(b_i - j),
    # where C(k) = 2 if k > 0 and C(0) = 1.
    # The sum is finite, running for j from 0 up to min(b).
    try:
        min_b = min(b)
    except ValueError:
        min_b = -1 # Handle empty list case

    total_invariant = 0
    contributions = []
    
    print("Calculating contributions for each term in the sum (from j=0 to j=min(b)):\n")
    for j in range(min_b + 1):
        exponents = [bi - j for bi in b]
        
        num_positive = sum(1 for p in exponents if p > 0)
        num_zero = n - num_positive # Since p cannot be negative in the sum
        
        # Contribution is 2^num_positive * 1^num_zero
        contribution = 2**num_positive
        contributions.append(contribution)
        total_invariant += contribution
        
        print(f"Term for j = {j}:")
        print(f"  Exponents (b_i - j): {exponents}")
        # The final equation requires each number. 2^{num_positive} shows the components.
        print(f"  Contribution = 2^{num_positive} * 1^{num_zero} = {contribution}")
        
    # Step 4: Print the final result and equation
    print("\n-----------------------------------------------------")
    print("The Crawley-Nordström invariant is the sum of these contributions.")
    
    equation_parts = [str(c) for c in contributions]
    equation_str = " + ".join(equation_parts)
    
    print(f"Final Equation: {equation_str} = {total_invariant}")
    print("-----------------------------------------------------")
    
    return total_invariant

# Run the calculation and store the final answer
result = calculate_crawley_nordstrom_invariant()
# The final answer will be enclosed in <<<>>> at the very end.