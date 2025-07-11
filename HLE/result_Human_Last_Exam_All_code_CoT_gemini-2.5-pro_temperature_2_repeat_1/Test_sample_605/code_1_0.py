import itertools

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau hypersurface
    in a weighted projective space.
    """
    # The weights of the ambient projective space P(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]
    
    # For a Calabi-Yau hypersurface, the degree 'd' of the defining polynomial
    # is the sum of the weights.
    d = sum(weights)
    
    # The invariant is given by the formula: d * sum(w_i * w_j for i < j).
    # We calculate the sum of pairwise products of the weights.
    sum_of_pairwise_products = 0
    for i, j in itertools.combinations(range(len(weights)), 2):
        sum_of_pairwise_products += weights[i] * weights[j]
        
    # The final invariant is the product of the degree and the sum of pairwise products.
    invariant = d * sum_of_pairwise_products
    
    # Output the final equation as requested.
    print(f"The weights are: {weights}")
    print(f"The degree of the Calabi-Yau hypersurface is d = sum(weights) = {d}")
    print(f"The sum of pairwise products of weights is sum(w_i*w_j) = {sum_of_pairwise_products}")
    print("\nThe Crawley-Nordstrom invariant is d * sum(w_i*w_j).")
    print("\nFinal Equation:")
    print(f"{d} * {sum_of_pairwise_products} = {invariant}")

calculate_crawley_nordstrom_invariant()