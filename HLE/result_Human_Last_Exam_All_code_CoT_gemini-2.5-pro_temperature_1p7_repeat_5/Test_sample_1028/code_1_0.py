def calculate_model_parameters():
    """
    Calculates the number of parameters for two different birth-death models
    to illustrate why hyper-complex models do not mitigate identifiability issues.
    """
    num_pieces = 10
    poly_degree = 5

    # --- Analysis for Strategy C: 10 pieces with polynomials of degree 5 ---
    # A polynomial of degree 'd' has d+1 coefficients (e.g., ax^2+bx+c has 3).
    coeffs_per_poly = poly_degree + 1
    
    # We need one set of coefficients for the speciation rate function
    # and another for the extinction rate function, for each time piece.
    params_for_C = num_pieces * coeffs_per_poly * 2  # x2 for speciation AND extinction

    print("--- Strategy C Analysis ---")
    print(f"Model: Birth-death with {num_pieces} time pieces.")
    print(f"Rates in each piece are defined by a polynomial of degree {poly_degree}.")
    print(f"Number of coefficients per polynomial: {poly_degree} (degree) + 1 = {coeffs_per_poly}")
    print(f"Number of rate functions per piece: 2 (1 for speciation, 1 for extinction)")
    print("Total rate parameters to estimate:")
    print(f"{num_pieces} (pieces) * {coeffs_per_poly} (coeffs per rate) * 2 (rates) = {params_for_C}")
    print("\nThis model has a massive number of parameters that cannot be identified from extant-only data, making the identifiability problem worse.\n")

    # --- Analysis for Strategy A: 10 pieces with constant rates ---
    # A constant rate is a polynomial of degree 0, so it has 1 coefficient.
    coeffs_for_constant = 0 + 1
    params_for_A = num_pieces * coeffs_for_constant * 2

    print("--- Strategy A Analysis (for comparison) ---")
    print(f"Model: Birth-death with {num_pieces} constant pieces.")
    print("Each constant rate is a polynomial of degree 0.")
    print("Number of coefficients per constant rate: 0 (degree) + 1 = 1")
    print(f"Number of rate functions per piece: 2 (1 for speciation, 1 for extinction)")
    print("Total rate parameters to estimate:")
    print(f"{num_pieces} (pieces) * {coeffs_for_constant} (coeffs per rate) * 2 (rates) = {params_for_A}")
    print("\nWhile this model also has unidentifiable parameters, it is far simpler than the model in C.")


calculate_model_parameters()