import math

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for the given Calabi-Yau Link.
    """
    # The weights of the ambient space
    weights = [22, 29, 49, 50, 75]

    # The polynomial is f = z_1^8*z_3 + ... + z_5^3.
    # For a polynomial to be quasi-homogeneous, all its monomials must have the same weighted degree.
    # Let's calculate the degree 'd' using the last monomial, z_5^3:
    # d = 3 * w_5 = 3 * 75 = 225.
    d = 225

    # The "Calabi-Yau" condition requires that the degree 'd' equals the sum of the weights.
    # Let's calculate the sum of the weights.
    sum_of_weights = sum(weights)

    # Let's calculate the product of the weights.
    product_of_weights = math.prod(weights)

    # The formula for the Crawley-Nordström invariant is CN = (d - sum(w_i)) / prod(w_i).
    # Since the Calabi-Yau condition (d = sum(w_i)) is met, the numerator will be 0.
    numerator = d - sum_of_weights
    
    # Calculate the invariant.
    if product_of_weights == 0:
        invariant = float('inf') if numerator != 0 else 0 # Avoid division by zero
    else:
        invariant = numerator / product_of_weights

    # Print the step-by-step calculation
    print("Step 1: Define the degree and weights.")
    print(f"Weights (w_i): {weights}")
    print(f"Degree (d): {d}")
    
    print("\nStep 2: Check the Calabi-Yau condition (d = sum(w_i)).")
    print(f"Sum of weights = {' + '.join(map(str, weights))} = {sum_of_weights}")
    print(f"Since d ({d}) = sum(w_i) ({sum_of_weights}), the Calabi-Yau condition is satisfied.")

    print("\nStep 3: Calculate the Crawley-Nordström invariant using the formula.")
    print("CN = (d - sum(w_i)) / prod(w_i)")
    print(f"Product of weights = {' * '.join(map(str, weights))} = {product_of_weights}")
    
    print("\nFinal Equation:")
    print(f"CN = ({d} - {sum_of_weights}) / {product_of_weights}")
    print(f"CN = {numerator} / {product_of_weights}")
    print(f"CN = {invariant}")

solve_crawley_nordstrom()
<<<0.0>>>