def solve_neron_severi_rank():
    """
    Calculates the smallest and largest possible ranks of the Neron-Severi group
    of the 15th symmetric power of a genus 3 Riemann surface.
    """

    # 1. Define the properties of the Riemann surface C and its symmetric power X.
    genus_g = 3
    degree_d = 15

    print(f"Let C be a Riemann surface of genus g = {genus_g}.")
    print(f"Let X be its {degree_d}th symmetric power, X = C^({degree_d}).")
    print("We want to find the range of the rank of the Neron-Severi group of X, denoted rho(X).")
    print("-" * 30)

    # 2. State the key formula.
    print("The rank rho(X) is related to the rank of the Neron-Severi group of the Jacobian of C, rho(J(C)), by the formula:")
    print(f"rho(C^({degree_d})) = 1 + rho(J(C))")
    print(f"This formula is valid because the degree d={degree_d} is greater than or equal to g-1 = {genus_g - 1}.")
    print("-" * 30)

    # 3. Determine the range for rho(J(C)).
    print("To find the range of rho(X), we need the range of rho(J(C)) for a genus 3 curve.")

    # Smallest rank (for a generic curve)
    rho_J_C_min = 1
    print(f"\nThe smallest possible rank, for a generic curve C, is rho(J(C)) = {rho_J_C_min}.")

    # Largest rank (for a curve with maximal Complex Multiplication)
    rho_J_C_max = genus_g**2
    print(f"The largest possible rank, for a special curve with maximal Complex Multiplication, is rho(J(C)) = g^2 = {genus_g}^2 = {rho_J_C_max}.")
    print("-" * 30)

    # 4. Calculate the smallest and largest ranks for X.
    print("Now we can calculate the smallest and largest possible ranks for X.")

    # Smallest rank
    smallest_rank = 1 + rho_J_C_min
    print(f"Smallest rank of NS(X) = 1 + rho(J(C))_min = 1 + {rho_J_C_min} = {smallest_rank}")

    # Largest rank
    largest_rank = 1 + rho_J_C_max
    print(f"Largest rank of NS(X)  = 1 + rho(J(C))_max = 1 + {rho_J_C_max} = {largest_rank}")

if __name__ == '__main__':
    solve_neron_severi_rank()