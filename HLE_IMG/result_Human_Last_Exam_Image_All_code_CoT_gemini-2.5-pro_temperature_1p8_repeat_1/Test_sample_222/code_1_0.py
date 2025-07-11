import math
from fractions import Fraction

def calculate_and_print_solution():
    """
    This function calculates the imaginary part of the sum of the contour integrals as described.
    It follows the step-by-step plan outlined above.
    """

    # Step 2: Define a function to calculate residues at poles z = -n
    def residue_at_neg_n(n):
        """Calculates the residue of f(z) at the pole z = -n."""
        # Res(f, -n) = (2*n / (2*n+3)) * ((-1)^n / n!)
        numerator = 2 * n * ((-1)**n)
        denominator = (2 * n + 3) * math.factorial(n)
        return Fraction(numerator, denominator)

    # Calculate residues for the poles with non-zero total winding numbers
    # From Step 5, we only need Res(f, -2) and Res(f, -3).
    res_minus_2 = residue_at_neg_n(2)
    res_minus_3 = residue_at_neg_n(3)
    
    # Step 5 & 6: Combine everything using the Residue Theorem
    # I = 2*pi*i * ( Ind(-2)*Res(-2) + Ind(-3)*Res(-3) )
    # Ind(-2) = 2, Ind(-3) = -1
    
    sum_of_weighted_residues = 2 * res_minus_2 - 1 * res_minus_3
    
    # The sum of integrals I is purely imaginary: I = i * (2 * pi * sum_of_weighted_residues)
    # The imaginary part is 2 * pi * sum_of_weighted_residues.

    # Print the explanation and the step-by-step calculation
    print("The total integral I is given by the Residue Theorem:")
    print("I = 2*pi*i * [ Ind_C(-2)*Res(f, -2) + Ind_C(-3)*Res(f, -3) ]")
    print("From the contours, we determined the total winding numbers: Ind_C(-2) = 2 and Ind_C(-3) = -1.\n")

    print("The residues are calculated as:")
    print(f"Res(f, -2) = 2*2/ (2*2+3) * (-1)^2/2! = 4/7 * 1/2 = {res_minus_2}")
    print(f"Res(f, -3) = 2*3/ (2*3+3) * (-1)^3/3! = 6/9 * -1/6 = {res_minus_3}\n")
    
    print("The sum of the weighted residues is:")
    print(f"Sum = 2 * Res(f, -2) - 1 * Res(f, -3)")
    print(f"    = 2 * ({res_minus_2}) - ({res_minus_3})")
    print(f"    = {2 * res_minus_2} + {(-1) * res_minus_3}")
    print(f"    = {sum_of_weighted_residues}\n")

    print("The imaginary part of the sum of the integrals is:")
    final_eq = f"Im(I) = 2 * pi * ( {sum_of_weighted_residues} )"
    print(final_eq)

    final_symbolic_val = f"{2 * sum_of_weighted_residues.numerator}/{sum_of_weighted_residues.denominator} * pi"
    final_numeric_val = float(2 * math.pi * sum_of_weighted_residues)
    print(f"Im(I) = {final_symbolic_val}")
    print(f"Im(I) approx. {final_numeric_val:.6f}")


calculate_and_print_solution()
<<<86/63 * pi>>>