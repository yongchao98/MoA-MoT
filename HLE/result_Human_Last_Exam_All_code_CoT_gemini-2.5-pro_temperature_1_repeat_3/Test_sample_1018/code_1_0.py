import math

def solve_simplicial_volume():
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.
    """
    # Step 1: Define the genera of the surfaces.
    g1 = 31
    g2 = 17

    # The dimensions of the surfaces are both 2.
    m = 2
    n = 2

    print(f"We want to compute the simplicial volume of Σ_{g1} x Σ_{g2}.")
    print("-" * 30)

    # Step 2: Calculate the simplicial volume of each surface.
    # The formula is ||Σ_g|| = 4 * (g - 1) for g >= 1.
    sv1 = 4 * (g1 - 1)
    sv2 = 4 * (g2 - 1)
    print(f"The simplicial volume of the first surface, ||Σ_{g1}||, is calculated as:")
    print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = 4 * {g1 - 1} = {sv1}")
    print()
    print(f"The simplicial volume of the second surface, ||Σ_{g2}||, is calculated as:")
    print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = 4 * {g2 - 1} = {sv2}")
    print("-" * 30)

    # Step 3: Calculate the binomial coefficient for the product formula.
    # The product formula is ||M x N|| = C(m+n, m) * ||M|| * ||N||.
    coefficient = math.comb(m + n, m)
    print("For the product space, we use Gromov's product formula:")
    print("||M x N|| = C(m+n, m) * ||M|| * ||N||")
    print(f"Here, the dimensions are m={m} and n={n}.")
    print(f"The coefficient is C({m}+{n}, {m}) = C({m+n}, {m}) = {coefficient}.")
    print("-" * 30)

    # Step 4: Calculate the final simplicial volume of the product space.
    total_sv = coefficient * sv1 * sv2
    print("Finally, we multiply the components to get the total simplicial volume:")
    print(f"||Σ_{g1} x Σ_{g2}|| = C({m+n}, {m}) * ||Σ_{g1}|| * ||Σ_{g2}||")
    print(f"The final equation with the computed values is:")
    print(f"{coefficient} * {sv1} * {sv2} = {total_sv}")

solve_simplicial_volume()