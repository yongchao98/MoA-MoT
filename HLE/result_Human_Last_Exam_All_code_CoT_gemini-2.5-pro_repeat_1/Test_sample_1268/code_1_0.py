import math

def solve_norm_covolume_bound():
    """
    Calculates the coefficient for the upper bound of the maximum norm 
    in relation to the covolume for quadratic fields and prints the resulting inequality.
    """
    
    # General parameters for quadratic fields
    n = 2  # Degree of the field

    # --- Case 1: Real Quadratic Field (e.g., Q(sqrt(N)), N > 0) ---
    # For a real quadratic field, there are 2 real embeddings and 0 complex embeddings.
    r2_real = 0
    
    # The Minkowski bound is Norm(a) <= ( (4/pi)^r2 * n! / n^n ) * sqrt(|D_K|)
    # Let's calculate the first part of the coefficient.
    minkowski_const_real = ((4 / math.pi)**r2_real * math.factorial(n)) / (n**n)
    
    # The covolume V is V = 2^(-r2) * sqrt(|D_K|).
    # So, sqrt(|D_K|) = 2^r2 * V.
    # The relationship is Norm(a) <= minkowski_const * (2^r2 * V).
    # The final coefficient C is minkowski_const * 2^r2.
    C_real = minkowski_const_real * (2**r2_real)

    # --- Case 2: Imaginary Quadratic Field (e.g., Q(sqrt(-N)), N > 0) ---
    # For an imaginary quadratic field, there are 0 real embeddings and 1 pair of complex embeddings.
    r2_imaginary = 1

    # Calculate the coefficient C for this case using the same logic.
    minkowski_const_imaginary = ((4 / math.pi)**r2_imaginary * math.factorial(n)) / (n**n)
    C_imaginary = minkowski_const_imaginary * (2**r2_imaginary)

    # The general upper bound is determined by the maximum of the two coefficients.
    # C_real = (1 * 2 / 4) * 1 = 0.5
    # C_imaginary = ((4/pi) * 2 / 4) * 2 = (2/pi) * 2 = 4/pi ~= 1.273
    # The larger coefficient comes from the imaginary case.
    
    # We will present the final answer symbolically as it is more precise.
    # The maximum coefficient is 4/π.
    norm_bound_symbol = "k_{k,∞}"
    covolume_symbol = "V"
    numerator = 4
    denominator_symbol = "π"

    print("The relationship between the maximum norm ({}) and the covolume ({}) for any quadratic field is given by the inequality:".format(norm_bound_symbol, covolume_symbol))
    
    # The final equation, with each number (4 and π) included.
    print()
    print("{} ≤ ({}/{}) * {}".format(norm_bound_symbol, numerator, denominator_symbol, covolume_symbol))
    print()
    
    # For verification, print the numeric values of the coefficients.
    print("This is because the coefficient for imaginary fields ({:.4f}) is larger than for real fields ({:.4f}).".format(C_imaginary, C_real))


solve_norm_covolume_bound()