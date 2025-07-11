import fractions

def calculate_energy_shift_coefficient(n, l):
    """
    Calculates the numerical coefficient for the second-order relativistic energy shift.

    The formula for the energy shift is:
    Delta_E = (E_R * alpha^4) * [3/(4*n^4) + 2/(n^3*(l+0.5)) - 3/(n^2*(l+0.5)^2)]
    where E_R = 1/2 * m * c^2 * alpha^2.
    This function calculates the bracketed term and the final coefficient for m*c^2*alpha^6.
    """
    
    # Use fractions for exact arithmetic
    n = fractions.Fraction(n)
    l = fractions.Fraction(l)
    
    l_plus_half = l + fractions.Fraction(1, 2)
    
    term1 = fractions.Fraction(3, 4 * n**4)
    term2 = fractions.Fraction(2, n**3 * l_plus_half)
    term3 = -fractions.Fraction(3, n**2 * l_plus_half**2)
    
    bracket_term = term1 + term2 + term3
    
    # Final coefficient for m*c^2*alpha^6
    # Delta_E = (1/2 * m * c^2 * alpha^2) * alpha^4 * bracket_term
    # Delta_E = (1/2 * bracket_term) * m * c^2 * alpha^6
    final_coeff = fractions.Fraction(1, 2) * bracket_term
    
    return final_coeff

# Given quantum numbers
n = 3
l = 2

# Calculate the coefficient
coeff = calculate_energy_shift_coefficient(n, l)

# Output the final result
# The final result is Delta_E = coeff * m * c^2 * alpha^6
numerator = coeff.numerator
denominator = coeff.denominator

print("The second-order energy shift is:")
if numerator > 0:
    print(f"Delta_E = ({numerator}/{denominator}) * m * c^2 * alpha^6")
else:
    # Print the negative sign separately for clarity as requested by the format
    print(f"Delta_E = -({abs(numerator)}/{denominator}) * m * c^2 * alpha^6")
    
# For the final answer block, let's print the parts of the equation
print("\nFinal equation parts:")
print(f"Sign: -")
print(f"Numerator: {abs(numerator)}")
print(f"Denominator: {denominator}")
print(f"Expression: m * c^2 * alpha^6")

final_answer = f"-({abs(numerator)}/{denominator}) * m * c^2 * alpha^6"
# <<<f"-({abs(coeff.numerator)}/{coeff.denominator}) * m * c^2 * alpha^6">>>