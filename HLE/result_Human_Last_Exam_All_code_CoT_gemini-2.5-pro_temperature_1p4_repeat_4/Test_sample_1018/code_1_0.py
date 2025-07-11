import math

def compute_simplicial_volume_product():
    """
    This function computes the simplicial volume of the product of two surfaces,
    Sigma_g1 and Sigma_g2, and prints the step-by-step calculation.
    """
    g1 = 31
    g2 = 17

    # Dimensions of the surfaces
    n1 = 2
    n2 = 2

    print(f"We want to compute the simplicial volume of Sigma_{g1} x Sigma_{g2}.")
    print("-----------------------------------------------------------------")

    # Step 1: Calculate the simplicial volume of the first surface
    # The formula for a surface of genus g >= 2 is ||Sigma_g|| = 4g - 4
    sv1 = 4 * g1 - 4
    print(f"First, we compute the simplicial volume of the surface Sigma_{g1}.")
    print(f"Using the formula ||Sigma_g|| = 4*g - 4 for g >= 2:")
    print(f"||Sigma_{g1}|| = 4 * {g1} - 4 = {sv1}")
    print("-----------------------------------------------------------------")

    # Step 2: Calculate the simplicial volume of the second surface
    sv2 = 4 * g2 - 4
    print(f"Next, we compute the simplicial volume of the surface Sigma_{g2}.")
    print(f"||Sigma_{g2}|| = 4 * {g2} - 4 = {sv2}")
    print("-----------------------------------------------------------------")

    # Step 3: Calculate the binomial coefficient for the product formula
    # The formula is ||M1 x M2|| = C(n1+n2, n1) * ||M1|| * ||M2||
    total_dim = n1 + n2
    coefficient = math.comb(total_dim, n1)
    print("For the product M1 x M2, the simplicial volume is given by:")
    print("||M1 x M2|| = C(dim(M1)+dim(M2), dim(M1)) * ||M1|| * ||M2||")
    print(f"Here, the dimensions are n1={n1} and n2={n2}.")
    print(f"The coefficient is C({n1}+{n2}, {n1}) = C({total_dim}, {n1}) = {coefficient}.")
    print("-----------------------------------------------------------------")

    # Step 4: Calculate the final result
    total_sv = coefficient * sv1 * sv2
    print("Finally, we multiply the three parts together:")
    print(f"||Sigma_{g1} x Sigma_{g2}|| = {coefficient} * ||Sigma_{g1}|| * ||Sigma_{g2}||")
    print(f"The final equation is: {coefficient} * {sv1} * {sv2} = {total_sv}")

compute_simplicial_volume_product()