def solve_surface_genus():
    """
    Calculates the smallest possible genus for a closed surface containing
    a given surface with genus 10 and one unknotted boundary component.
    """
    # The given genus of the initial surface Sigma.
    g_Sigma = 10

    print(f"The genus of the initial surface Σ is g(Σ) = {g_Sigma}.")
    print("The final closed surface Σ' is formed by capping the boundary of Σ with a surface D.")
    print("The genus of the resulting surface follows the formula: g(Σ') = g(Σ) + g(D).\n")

    # To find the smallest possible genus for Σ', we need the smallest possible genus for D.
    # The genus of a surface must be a non-negative integer.
    # The smallest possible genus is 0.
    g_D_min = 0

    print(f"To minimize g(Σ'), we must minimize g(D).")
    print(f"The minimum possible genus for any surface is {g_D_min}.")
    print("A surface D with genus 0 and one boundary is a disk.")
    print("Since the boundary of Σ is unknotted, it is guaranteed that it can be capped by an embedded disk.\n")
    print("Therefore, the minimum genus for the capping surface D is g(D)_min = 0.")

    # Calculate the minimum possible genus for the final surface Σ'.
    g_Sigma_prime_min = g_Sigma + g_D_min

    # Print the final calculation.
    print("\nCalculating the minimum genus for the final surface Σ':")
    print(f"g(Σ')_min = g(Σ) + g(D)_min")
    print(f"g(Σ')_min = {g_Sigma} + {g_D_min}")
    print(f"g(Σ')_min = {g_Sigma_prime_min}")

    print(f"\nThis result holds regardless of the specific choice of Σ, as long as it meets the given criteria.")
    print(f"Any other choice for the capping surface D would have a genus g(D) > 0, leading to a larger g(Σ').")
    print(f"Thus, the smallest positive integer g is {g_Sigma_prime_min}.")

solve_surface_genus()