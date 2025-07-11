import math

def compute_simplicial_volume_product():
    """
    Computes the simplicial volume of the product of two surfaces, Sigma_31 and Sigma_17.
    """
    # Step 0: Define the genera and dimensions
    g1 = 31
    g2 = 17
    dim1 = 2
    dim2 = 2

    print(f"We want to compute the simplicial volume of Sigma_{g1} x Sigma_{g2}.")
    print("The formula is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||.")
    print("-" * 20)

    # Step 1: Compute the simplicial volume of the first surface, Sigma_31
    # The formula for a surface with genus g >= 2 is ||Sigma_g|| = 4g - 4.
    if g1 >= 2:
        sv1 = 4 * g1 - 4
        print(f"First, we compute the simplicial volume of Sigma_{g1}:")
        print(f"||Sigma_{g1}|| = 4 * {g1} - 4 = {sv1}")
    else:
        sv1 = 0
        print(f"The simplicial volume of Sigma_{g1} is {sv1} because its genus is less than 2.")
    print("-" * 20)

    # Step 2: Compute the simplicial volume of the second surface, Sigma_17
    if g2 >= 2:
        sv2 = 4 * g2 - 4
        print(f"Next, we compute the simplicial volume of Sigma_{g2}:")
        print(f"||Sigma_{g2}|| = 4 * {g2} - 4 = {sv2}")
    else:
        sv2 = 0
        print(f"The simplicial volume of Sigma_{g2} is {sv2} because its genus is less than 2.")
    print("-" * 20)

    # Step 3: Compute the binomial coefficient
    total_dim = dim1 + dim2
    C = math.comb(total_dim, dim1)
    print(f"Then, we compute the binomial coefficient C({total_dim}, {dim1}):")
    print(f"C({dim1} + {dim2}, {dim1}) = {C}")
    print("-" * 20)

    # Step 4: Compute the final result
    final_volume = C * sv1 * sv2
    print("Finally, we multiply these values together:")
    print(f"||Sigma_{g1} x Sigma_{g2}|| = {C} * {sv1} * {sv2} = {final_volume}")

compute_simplicial_volume_product()