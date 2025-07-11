import math

def solve_mass_ratio_cpn1():
    """
    This function determines the leading-order asymptotic mass ratio between the
    lightest and the subsequent higher solitonic excitations in the spectrum of
    CP(N-1) models as N approaches infinity.

    The derivation is presented step-by-step.
    """
    print("This program calculates the required mass ratio based on the semi-classical analysis of the CP(N-1) model.")
    print("-" * 80)

    # Step 1: State the mass formula
    print("Step 1: The mass spectrum formula")
    print("From semi-classical quantization in the CP(N-1) model, the mass M_k of the k-th stable excitation is given by:")
    print("\n  M_k = C * sin(pi * k / N)\n")
    print("where C is a constant that depends on N and the dynamically generated mass scale, but not on the integer k (where k = 1, 2, ..., N-1).\n")

    # Step 2: Define the specific masses
    print("Step 2: Define the masses for the lightest (k=1) and next (k=2) excitations")
    print("The lightest solitonic excitation corresponds to k=1:")
    print("  M_1 = C * sin(pi * 1 / N)")
    print("\nThe subsequent higher excitation corresponds to k=2:")
    print("  M_2 = C * sin(pi * 2 / N)\n")

    # Step 3: Form the ratio
    print("Step 3: Form the mass ratio R = M_2 / M_1")
    print("The ratio R is independent of the constant C:")
    print("  R = M_2 / M_1 = [C * sin(2*pi/N)] / [C * sin(pi/N)]")
    print("  R = sin(2*pi/N) / sin(pi/N)\n")

    # Step 4: Simplify and evaluate the limit
    print("Step 4: Simplify the ratio and evaluate the limit as N -> infinity")
    print("Using the trigonometric identity sin(2*theta) = 2 * sin(theta) * cos(theta), with theta = pi/N, we get:")
    print("  R = [2 * sin(pi/N) * cos(pi/N)] / sin(pi/N)")
    print("After cancelling sin(pi/N) from the numerator and denominator, the expression simplifies to:")
    print("  R = 2 * cos(pi/N)\n")

    # Step 5: Final calculation
    print("Step 5: Calculate the final asymptotic value")
    print("To find the asymptotic ratio, we take the limit as N approaches infinity:")
    print("  lim (N->inf) R = lim (N->inf) [2 * cos(pi/N)]")
    print("As N -> inf, the term (pi/N) -> 0. Since cos(0) = 1, the limit is:")
    
    # Final equation with numbers
    final_ratio_coefficient = 2
    cos_limit_value = 1
    final_result = final_ratio_coefficient * cos_limit_value
    print(f"\nFinal Asymptotic Ratio = {final_ratio_coefficient} * {cos_limit_value} = {final_result}")
    print("-" * 80)

if __name__ == "__main__":
    solve_mass_ratio_cpn1()

<<<2>>>