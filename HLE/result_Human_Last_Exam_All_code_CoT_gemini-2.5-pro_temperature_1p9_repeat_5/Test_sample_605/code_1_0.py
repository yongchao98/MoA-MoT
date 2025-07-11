def calculate_invariant():
    """
    Calculates the Crawley-Nordström invariant, interpreted as the Milnor number,
    for the given Calabi-Yau Link.
    """
    weights = [22, 29, 49, 50, 75]
    degree = sum(weights)
    
    # Check if the sum of weights is indeed 225
    # print(f"Weights: {weights}")
    # print(f"Degree d = sum of weights = {degree}")

    factors = []
    numerators = []
    denominators = []
    
    for w in weights:
        # Calculate each term (d/w - 1) which is (d-w)/w
        num = degree - w
        den = w
        factors.append(num / den)
        numerators.append(num)
        denominators.append(den)
    
    # Calculate the product of the factors for the final result
    result = 1
    for f in factors:
        result *= f
    
    # The result should be an integer, so we round it to avoid float inaccuracies
    final_result = round(result)

    # Building the equation string
    # (d/w_1-1) * (d/w_2-1) * ...
    # which is (d-w_1)/w_1 * (d-w_2)/w_2 * ...
    equation_parts = []
    for i in range(len(weights)):
        equation_parts.append(f"({numerators[i]}/{denominators[i]})")
        
    equation_str = " * ".join(equation_parts)
    
    print(f"The invariant is calculated by the formula μ = ∏(d/w_i - 1):")
    print(f"μ = {equation_str}")

    term_values = [f"{n/d:.4f}" for n,d in zip(numerators, denominators)]
    term_values_str = " * ".join(term_values)

    print(f"μ = {term_values_str}")
    print(f"μ = {final_result}")

calculate_invariant()