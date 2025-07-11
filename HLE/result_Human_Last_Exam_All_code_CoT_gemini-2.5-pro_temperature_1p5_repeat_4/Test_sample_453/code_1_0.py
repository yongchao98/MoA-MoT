import sympy

def calculate_relativistic_shift():
    """
    Calculates the first-order energy shift for the relativistic kinetic
    energy correction in a hydrogen atom.
    The calculation is performed symbolically using sympy.
    """
    # Define symbols for the fundamental constants and quantum numbers
    m, c, alpha, n, l = sympy.symbols('m c alpha n l', positive=True, real=True)
    
    print("Step 1: Define the first-order energy shift formula.")
    # The first-order shift is the expectation value of the perturbation H'
    # H' = -p^4 / (8 * m^3 * c^2)
    # We can write p^2 = 2m * T = 2m * (H0 - V)
    # H' = -(2m * (H0 - V))^2 / (8 * m^3 * c^2) = -1/(2*m*c**2) * (H0 - V)**2
    # The expectation value for state |n,l> is:
    # E_shift = <n,l| H' |n,l> = -1/(2*m*c**2) * <(E_n - V)^2>
    # E_shift = -1/(2*m*c**2) * (E_n**2 - 2*E_n*<V> + <V**2>)
    print("E_shift = -1/(2*m*c**2) * (E_n**2 - 2*E_n*<V>_nl + <V**2>_nl)\n")

    print("Step 2: Define the expressions for the terms.")
    # Unperturbed energy of the hydrogen atom
    E_n = -m * c**2 * alpha**2 / (2 * n**2)
    print(f"Unperturbed energy E_n = {E_n}")
    
    # Expectation value of the potential V, from the Virial Theorem <V> = 2*E_n
    V_exp = 2 * E_n
    print(f"Expectation value <V>_nl = {V_exp}")

    # Expectation value of V^2. V = -k_e*e^2/r = -alpha*hbar*c/r
    # <V^2> = (alpha*hbar*c)^2 * <1/r^2>
    # <1/r^2>_nl = 1 / (a0**2 * n**3 * (l + 1/2)), with a0 = hbar/(m*c*alpha)
    # So, <1/r^2>_nl = (m*c*alpha/hbar)**2 / (n**3 * (l + 1/2))
    # <V^2>_nl = (alpha*hbar*c)^2 * (m*c*alpha/hbar)**2 / (n**3 * (l + 1/2))
    # <V^2>_nl = m**2 * c**4 * alpha**4 / (n**3 * (l + 1/2))
    V2_exp = (m**2 * c**4 * alpha**4) / (n**3 * (l + sympy.S(1)/2))
    print(f"Expectation value <V^2>_nl = {V2_exp}\n")

    print("Step 3: Substitute these into the energy shift formula.")
    # E_shift = -1/(2*m*c**2) * (E_n**2 - 2*E_n*(2*E_n) + V2_exp)
    # E_shift = -1/(2*m*c**2) * (-3*E_n**2 + V2_exp)
    E_shift_expr = -1/(2*m*c**2) * (-3*E_n**2 + V2_exp)
    print("Substituting E_n and <V^2> into the formula for E_shift gives:")
    print(f"E_shift = {E_shift_expr}\n")
    
    print("Step 4: Simplify the expression.")
    E_shift_simplified = sympy.simplify(E_shift_expr)
    print(f"Simplified general formula for the shift:\nE_shift = {E_shift_simplified}\n")
    
    print("Step 5: Substitute n=3 and l=2 into the formula.")
    # Given quantum numbers are n=3, l=2
    final_shift = E_shift_simplified.subs({n: 3, l: 2})
    
    # Extract the numerator and denominator for clean printing
    num, den = final_shift.as_numer_denom()
    
    # We want to show the numbers in the equation m*c**2*alpha**4 * (coeff)
    # The non-constant part is m*c**2*alpha**4
    coeff = final_shift / (m * c**2 * alpha**4)
    coeff_num, coeff_den = coeff.as_numer_denom()
    
    print("For n=3, l=2, the energy shift is:")
    print(f"E_shift = ({coeff_num} / {coeff_den}) * m * c**2 * alpha**4")
    
    # Print each number in the final equation
    print(f"The numbers in the final equation are {abs(coeff_num)}, {coeff_den}.")

    return final_shift

# Run the calculation
final_answer = calculate_relativistic_shift()
# Format the final answer for the checker
final_answer_str = str(final_answer).replace('**', '^')
# We need to express it as - (fraction) * m * c^2 * alpha^4
coeff_final = final_answer / (sympy.symbols('m') * sympy.symbols('c')**2 * sympy.symbols('alpha')**4)
final_expression = f"{coeff_final}*m*c^2*alpha^4"

# Although the final answer is an expression, we are instructed to give a content answer
# The coefficient is -1/360
# Let's use the final symbolic expression as the answer content
# Final symbolic expression is -m*c**2*alpha**4/360
final_symbolic_expr_str = "-m*c**2*alpha**4/360"
print(f"\nFinal symbolic expression: {final_symbolic_expr_str}")

# Final Answer Block
# <<<answer content>>> format requires a single value or expression.
final_answer_for_checker = "-m*c**2*alpha**4/360"
print(f"\n<<<{final_answer_for_checker}>>>")
