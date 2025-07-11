import math

def solve_simplicial_volume():
    """
    Computes the simplicial volume of the product of two surfaces Σ_g1 and Σ_g2.
    """
    # Define the genera of the two surfaces
    g1 = 31
    g2 = 17

    # The dimension of a surface is 2
    n1 = 2
    n2 = 2

    # For a surface Σ_g, the simplicial volume is non-zero only if g >= 2.
    # If either g < 2, the simplicial volume of the product is 0.
    if g1 < 2 or g2 < 2:
        print("One of the surfaces has genus less than 2, so its simplicial volume is 0.")
        print("Therefore, the simplicial volume of the product is 0.")
        return

    # Calculate the simplicial volume of each surface using the formula ||Σ_g|| = 4 * (g - 1)
    sv1 = 4 * (g1 - 1)
    sv2 = 4 * (g2 - 1)

    # Calculate the binomial coefficient C(n1+n2, n1) for the product formula
    n_sum = n1 + n2
    try:
        # math.comb is available in Python 3.8+
        C = math.comb(n_sum, n1)
    except AttributeError:
        # Fallback for older Python versions
        C = math.factorial(n_sum) // (math.factorial(n1) * math.factorial(n2))

    # Calculate the total simplicial volume of the product
    total_sv = C * sv1 * sv2

    # Print the step-by-step calculation
    print("The simplicial volume of the product of surfaces Σ_g1 x Σ_g2 is computed as follows:")
    print("Formula: ||Σ_g1 x Σ_g2|| = C(n1+n2, n1) * ||Σ_g1|| * ||Σ_g2||, where n1=dim(Σ_g1), n2=dim(Σ_g2)\n")

    print("Step 1: Compute the simplicial volume of each surface.")
    print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = {sv1}")
    print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = {sv2}\n")

    print("Step 2: Compute the binomial coefficient.")
    print(f"C({n1}+{n2}, {n1}) = C({n_sum}, {n1}) = {C}\n")

    print("Step 3: Compute the final simplicial volume of the product.")
    print(f"||Σ_{g1} x Σ_{g2}|| = {C} * {sv1} * {sv2} = {total_sv}")

solve_simplicial_volume()