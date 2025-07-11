import math

def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces, Σ_g1 and Σ_g2.
    """
    # Step 1: Define the genera of the surfaces.
    g1 = 31
    g2 = 17

    print(f"We are computing the simplicial volume of Σ_{g1} x Σ_{g2}.")
    print("-" * 20)

    # Step 2: Calculate the simplicial volume of the first surface.
    # The formula for the simplicial volume of a hyperbolic surface (g >= 2) is ||Σ_g|| = 4g - 4.
    if g1 >= 2:
        sv1 = 4 * g1 - 4
        print(f"The simplicial volume of Σ_{g1} is ||Σ_{g1}|| = 4 * {g1} - 4 = {sv1}.")
    else:
        # This case is not applicable here but included for completeness.
        sv1 = 0
        print(f"The simplicial volume of Σ_{g1} is 0 since g < 2.")
    
    # Step 3: Calculate the simplicial volume of the second surface.
    if g2 >= 2:
        sv2 = 4 * g2 - 4
        print(f"The simplicial volume of Σ_{g2} is ||Σ_{g2}|| = 4 * {g2} - 4 = {sv2}.")
    else:
        sv2 = 0
        print(f"The simplicial volume of Σ_{g2} is 0 since g < 2.")
    
    print("-" * 20)
    
    # Step 4: Use the product formula for simplicial volume.
    # For aspherical manifolds M (dim m) and N (dim n), ||M x N|| = C(m+n, m) * ||M|| * ||N||.
    # For surfaces, the dimension is 2. So m=2, n=2.
    m = 2
    n = 2
    
    # Calculate the binomial coefficient C(m+n, m)
    binomial_coeff = math.comb(m + n, m)
    print(f"The product formula requires the binomial coefficient C({m}+{n}, {m}) = C({m+n}, {m}) = {binomial_coeff}.")
    print("-" * 20)
    
    # Step 5: Compute the final result.
    final_result = binomial_coeff * sv1 * sv2
    
    print("The final simplicial volume is calculated as follows:")
    # Print the final equation with all numbers
    print(f"||Σ_{g1} x Σ_{g2}|| = {binomial_coeff} * {sv1} * {sv2} = {final_result}")

compute_simplicial_volume_product()