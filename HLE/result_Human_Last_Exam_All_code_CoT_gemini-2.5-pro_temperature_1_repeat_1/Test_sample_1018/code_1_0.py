import math

def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces Σ_g1 and Σ_g2.
    """
    # Define the genera of the two surfaces
    g1 = 31
    g2 = 17

    # The dimension of a surface is 2
    dim = 2

    print(f"This script computes the simplicial volume of Σ_{g1} x Σ_{g2}.")
    print("-" * 50)

    # Step 1: Calculate the simplicial volume of each surface
    print("Step 1: Calculate the simplicial volume of each surface.")
    print("The formula for the simplicial volume of a surface of genus g >= 2 is ||Σ_g|| = 4 * (g - 1).\n")

    # For the first surface
    if g1 >= 2:
        sv1 = 4 * (g1 - 1)
        print(f"For the first surface, Σ_{g1}:")
        print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = 4 * {g1 - 1} = {sv1}\n")
    else:
        sv1 = 0
        print(f"The genus of the first surface Σ_{g1} is less than 2, so its simplicial volume is {sv1}.\n")

    # For the second surface
    if g2 >= 2:
        sv2 = 4 * (g2 - 1)
        print(f"For the second surface, Σ_{g2}:")
        print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = 4 * {g2 - 1} = {sv2}\n")
    else:
        sv2 = 0
        print(f"The genus of the second surface Σ_{g2} is less than 2, so its simplicial volume is {sv2}.\n")

    # Step 2: Apply the product formula
    print("Step 2: Apply the product formula for simplicial volumes.")
    print("The formula is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||.\n")

    # Calculate the binomial coefficient
    total_dim = dim + dim
    binom_coeff = math.comb(total_dim, dim)
    print(f"The dimensions are m = {dim} and n = {dim}.")
    print(f"The binomial coefficient is C({total_dim}, {dim}) = {binom_coeff}.\n")

    # Step 3: Compute the final result
    print("Step 3: Compute the final simplicial volume.")
    total_sv = binom_coeff * sv1 * sv2
    
    print("The simplicial volume of the product Σ_{g1} x Σ_{g2} is:")
    print(f"||Σ_{g1} x Σ_{g2}|| = {binom_coeff} * ||Σ_{g1}|| * ||Σ_{g2}||")
    print(f"               = {binom_coeff} * {sv1} * {sv2}")
    print(f"               = {total_sv}")
    print("-" * 50)
    print(f"The final result is: {total_sv}")

if __name__ == "__main__":
    compute_simplicial_volume_product()