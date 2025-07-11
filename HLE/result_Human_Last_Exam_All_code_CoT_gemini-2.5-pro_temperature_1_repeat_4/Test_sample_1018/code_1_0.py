import math

def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.
    """
    # Step 1: Define the genera of the surfaces
    g1 = 31
    g2 = 17
    
    # The dimensions of the surfaces are n=2 and m=2
    n = 2
    m = 2

    print(f"We want to compute the simplicial volume of Σ_{g1} × Σ_{g2}.")
    
    # Step 2: Calculate the simplicial volume of each surface
    # The formula for the simplicial volume of a surface Σ_g with g >= 2 is ||Σ_g|| = 4g - 4.
    sv1 = 4 * g1 - 4
    sv2 = 4 * g2 - 4
    
    print(f"First, we compute the simplicial volume of each surface individually.")
    print(f"For Σ_{g1}, the simplicial volume is ||Σ_{g1}|| = 4 * {g1} - 4 = {sv1}.")
    print(f"For Σ_{g2}, the simplicial volume is ||Σ_{g2}|| = 4 * {g2} - 4 = {sv2}.")
    
    # Step 3: Use the product formula for simplicial volumes
    # The formula is ||M x N|| = (1 / C(n+m, n)) * ||M|| * ||N||
    # where C is the binomial coefficient.
    binom_coeff_val = math.comb(n + m, n)
    
    print("\nNext, we use Gromov's product formula for simplicial volumes:")
    print(f"||Σ_{g1} × Σ_{g2}|| = (1 / C({n+m}, {n})) * ||Σ_{g1}|| * ||Σ_{g2}||")
    print(f"The binomial coefficient C({n+m}, {n}) is {binom_coeff_val}.")
    
    # Step 4: Compute the final value
    total_sv = (sv1 * sv2) / binom_coeff_val
    
    print("\nPutting it all together, the final calculation is:")
    print(f"||Σ_{g1} × Σ_{g2}|| = (1 / {binom_coeff_val}) * {sv1} * {sv2} = {int(total_sv)}")

compute_simplicial_volume_product()