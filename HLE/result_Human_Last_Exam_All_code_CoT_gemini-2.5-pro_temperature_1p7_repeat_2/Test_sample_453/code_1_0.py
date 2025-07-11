import sympy

def calculate_relativistic_shift_second_order(n, l):
    """
    Calculates the second-order relativistic kinetic energy shift for a hydrogenic atom.

    Args:
        n: Principal quantum number.
        l: Angular momentum quantum number.
    """
    # Use sympy for precise fractional arithmetic
    n, l = sympy.S(n), sympy.S(l)

    if l == 0:
        # The formula for l=0 is different, but the problem specifies l=2.
        # This implementation is for l > 0 as per the source formula.
        print("This formula is for l > 0.")
        return

    # Calculate each term inside the bracket of the formula
    term1_base = 4 * n / (l + sympy.S(1)/2) - 3
    term1 = -(term1_base**2)
    
    term2 = -2 / (2 * l + 1)
    
    term3_num = 4 * l * (l + 1)
    term3_den = (2 * l - 1) * (2 * l + 1) * (2 * l + 3)
    term3 = term3_num / term3_den

    bracket_value = term1 + term2 + term3

    # Now, let's print the step-by-step calculation of the bracket
    print("Step 1: Calculate the terms inside the square brackets.")
    print(f"Given n = {n}, l = {l}:")
    print(f"Term 1: -((4*n / (l + 1/2)) - 3)^2 = -((4*{n} / ({l} + 1/2)) - 3)^2 = -({term1_base})^2 = {term1}")
    print(f"Term 2: -2 / (2*l + 1) = -2 / (2*{l} + 1) = {term2}")
    print(f"Term 3: (4*l*(l+1)) / ((2*l-1)*(2*l+1)*(2*l+3)) = {term3_num}/{term3_den} = {term3}")
    print(f"Sum of terms in bracket: {term1} + {term2} + {term3} = {bracket_value}")
    
    # Calculate the prefactor
    # The full prefactor is (m*c**2 * alpha**6) / (4 * n**4)
    prefactor_fraction = sympy.S(1) / (4 * n**4)
    
    final_fraction = prefactor_fraction * bracket_value
    
    print("\nStep 2: Combine with the prefactor (m*c^2 * alpha^6) / (4 * n^4).")
    print(f"Prefactor fraction: 1 / (4 * {n}^4) = {prefactor_fraction}")
    print(f"Final Energy Shift Delta_E(2) = ({prefactor_fraction}) * ({bracket_value}) * m*c^2*alpha^6")
    print("\nFinal Result:")
    print(f"Delta_E(2) = {final_fraction} * m * c^2 * alpha^6")

    # To show the numbers in the final equation as requested.
    num, den = final_fraction.as_numer_denom()
    print(f"\nWriting out the full expression for the energy shift:")
    print(f"Delta_E(2) = ({num} / {den}) * m * c^2 * alpha^6")


# Given values from the problem
n_val = 3
l_val = 2

calculate_relativistic_shift_second_order(n_val, l_val)

# For the final answer extraction
n_s, l_s = sympy.S(n_val), sympy.S(l_val)
term1_s = -((4 * n_s / (l_s + sympy.S(1)/2) - 3)**2)
term2_s = -2 / (2 * l_s + 1)
term3_s = (4 * l_s * (l_s + 1)) / ((2 * l_s - 1) * (2 * l_s + 1) * (2 * l_s + 3))
bracket_value_s = term1_s + term2_s + term3_s
prefactor_fraction_s = sympy.S(1) / (4 * n_s**4)
final_fraction_s = prefactor_fraction_s * bracket_value_s
final_num, final_den = final_fraction_s.as_numer_denom()
answer_str = f"({final_num}/{final_den}) * m_e * c^2 * alpha^6"
print(f'<<<{answer_str}>>>')