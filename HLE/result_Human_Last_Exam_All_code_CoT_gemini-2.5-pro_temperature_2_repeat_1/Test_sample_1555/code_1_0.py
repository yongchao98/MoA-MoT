import math

def solve_mass_ratio():
    """
    Determines the leading-order asymptotic mass ratio between the lightest
    and the next-to-lightest excitation in CP(N-1) models.
    """
    print("This script determines the asymptotic mass ratio in CP(N-1) models in the large-N limit.")
    print("-" * 70)

    # Step 1: Explain the physics and the mass formula.
    print("1. The Model and Mass Spectrum")
    print("In 2D CP(N-1) models, the mass spectrum of stable particles is given by a formula derived from the model's exact S-matrix.")
    print("The mass 'm_k' of the k-th bound state (for k=1, 2, ...) is given by:")
    print("\n  m_k = C * sin(k * pi / N)\n")
    print("Where C is a normalization constant, and N is the rank of the model.")
    print("-" * 70)

    # Step 2: Formulate and simplify the mass ratio.
    print("2. The Mass Ratio Calculation")
    print("We need the ratio between the lightest excitation (k=1) and the subsequent higher excitation (k=2).")
    print("  Ratio = m_2 / m_1")
    print("  Ratio = [C * sin(2 * pi / N)] / [C * sin(1 * pi / N)]")
    print("  Ratio = sin(2 * pi / N) / sin(pi / N)")
    print("\nUsing the trigonometric identity sin(2x) = 2*sin(x)*cos(x), the expression simplifies to:")
    print("  Ratio = [2 * sin(pi / N) * cos(pi / N)] / sin(pi / N)")
    print("  Ratio = 2 * cos(pi / N)")
    print("-" * 70)

    # Step 3: Evaluate the asymptotic limit.
    print("3. Asymptotic Limit as N approaches Infinity")
    print("We evaluate the limit of the ratio as N -> infinity:")
    print("  lim (N->inf) [2 * cos(pi / N)]")

    # The numbers in the final limiting equation are 2 and the limit of cos(pi/N).
    constant_factor = 2
    # As N -> infinity, (pi / N) -> 0. Therefore, cos(pi/N) -> cos(0).
    limit_of_cos_term = math.cos(0)

    final_ratio = constant_factor * limit_of_cos_term

    print(f"\nAs N approaches infinity, the term (pi / N) approaches 0.")
    print(f"The limit of cos(pi / N) is cos(0), which evaluates to: {limit_of_cos_term}")
    print("\nThe final equation for the limit is the product of these numerical values:")
    
    # Printing each number in the final equation as requested.
    print(f"  Final Equation: {constant_factor} * {limit_of_cos_term} = {final_ratio}")
    print("-" * 70)

    print(f"\nTherefore, the leading-order asymptotic mass ratio is {final_ratio}.")

if __name__ == "__main__":
    solve_mass_ratio()
<<<2>>>