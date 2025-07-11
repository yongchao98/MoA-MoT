import math

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau threefold
    hypersurface in a weighted projective space.

    The formula used is cn(X) = (1/12) * (sigma_2 * d / prod_w),
    where:
    - w are the weights of the ambient space
    - d is the degree of the hypersurface (sum of weights)
    - sigma_2 is the second elementary symmetric polynomial of the weights
    - prod_w is the product of the weights
    """
    # The weights of the ambient space P^4_(w1, w2, w3, w4, w5)
    weights = [22, 29, 49, 50, 75]

    # The degree 'd' for a Calabi-Yau is the sum of the weights
    d = sum(weights)

    # The product of the weights
    prod_w = math.prod(weights)

    # Calculate sigma_2 = sum_{i<j} w_i * w_j
    # We use the identity: 2 * sigma_2 = (sum(w_i))^2 - sum(w_i^2)
    sum_sq_w = sum(w**2 for w in weights)
    sigma_2 = (d**2 - sum_sq_w) / 2
    
    # Although sigma_2 must be an integer, ensure it is for the calculation
    sigma_2 = int(sigma_2)

    # Calculate the value of the integral of c_2 . H
    c2_dot_L = sigma_2 * (d / prod_w)
    
    # The Crawley-Nordstrom invariant
    invariant = c2_dot_L / 12

    # Output the components of the calculation and the final result
    print("The Calabi-Yau threefold is defined by weights w = (22, 29, 49, 50, 75).")
    print("The Crawley-Nordström invariant is calculated using the formula: cn(X) = (1/12) * sigma_2 * (d / product(w_i))")
    print("\nComponent values:")
    print(f"Degree d = sum(w_i) = {d}")
    print(f"sigma_2 = sum(w_i*w_j for i<j) = {sigma_2}")
    print(f"Product of weights = product(w_i) = {prod_w}")

    print("\nFinal Equation:")
    print(f"cn(X) = (1 / 12) * {sigma_2} * ({d} / {prod_w})")
    print(f"cn(X) = {invariant}")


if __name__ == '__main__':
    calculate_crawley_nordstrom_invariant()
    # The format below is for the final answer extraction.
    weights = [22, 29, 49, 50, 75]
    d = sum(weights)
    prod_w = math.prod(weights)
    sum_sq_w = sum(w**2 for w in weights)
    sigma_2 = int((d**2 - sum_sq_w) / 2)
    c2_dot_L = sigma_2 * (d / prod_w)
    invariant = c2_dot_L / 12
    print(f"\n<<<{invariant}>>>")