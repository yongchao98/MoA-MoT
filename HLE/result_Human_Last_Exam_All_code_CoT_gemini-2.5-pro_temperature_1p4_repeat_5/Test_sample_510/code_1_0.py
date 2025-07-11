import sympy

def solve_and_format():
    """
    This function calculates the optimal mean response time 'x' for the described
    queueing system, filters the expression based on the given rules, and prints
    the remaining term in the specified LaTeX format.
    """
    # Step 1: Define symbols and parameters
    s, y = sympy.symbols('s y', real=True, positive=True)
    lambda_val = sympy.Rational(3, 2)
    # For U(0,1), the PDF f(s) = 1 over the interval [0, 1]
    f_s = 1

    # Step 2: Calculate the mean service time E[S]
    E_S = sympy.integrate(s * f_s, (s, 0, 1))

    # Step 3: Define the components for the mean waiting time E[W] formula
    # ρ_s = λ * ∫[0,s] y*f(y) dy
    rho_s_integral = sympy.integrate(y * f_s, (y, 0, s))
    rho_s = lambda_val * rho_s_integral

    # E[S_s^2] = ∫[0,s] y^2*f(y) dy (This is the partial second moment)
    moment2_s_integral = sympy.integrate(y**2 * f_s, (y, 0, s))

    # Integrand for the E[W] calculation, where the outer integral is over s from 0 to 1
    integrand = (lambda_val * moment2_s_integral) / (2 * (1 - rho_s)**2) * f_s

    # Step 4: Calculate the mean waiting time E[W]
    E_W = sympy.integrate(integrand, (s, 0, 1))

    # Step 5: Calculate the optimal mean response time x = E[T] = E[S] + E[W]
    x_optimal = E_S + E_W

    # Step 6: Filter the expression for x
    # Remove rational terms and terms that are logarithms of rational numbers.
    remaining_term_sum = 0
    # The .as_ordered_terms() method splits the sum into its components
    terms = x_optimal.as_ordered_terms()

    for term in terms:
        # Check if the term is a rational number
        is_rational_term = term.is_Rational

        # Check if the term is a logarithm of a rational number, i.e., of the form log(q)
        is_log_of_rational_term = False
        if term.func == sympy.log:
            if len(term.args) == 1 and term.args[0].is_Rational:
                is_log_of_rational_term = True

        # Keep the term if it doesn't fall into the categories to be removed
        if not is_rational_term and not is_log_of_rational_term:
            remaining_term_sum += term

    # Step 7: Format the final answer in LaTeX
    # The remaining term is expected to be of the form: (rational_coeff) * log(rational_arg)
    final_term = remaining_term_sum
    
    print("The final equation for the remaining term contains the following numbers:")
    # Extract numbers for the "output each number" requirement.
    # We parse the remaining term which is -4/9*log(2)
    if final_term != 0:
        coeff, log_part = final_term.as_coeff_Mul()
        num, den = coeff.p, coeff.q
        log_arg = log_part.args[0]
        # Numbers are absolute value of numerator, denominator, and log argument.
        print(abs(num))
        print(den)
        print(log_arg)
        
        # Construct the LaTeX string as per the specified format
        latex_string = f"{'- ' if num < 0 else ''}\\frac{{{abs(num)}}}{{{den}}} \\ln({{{log_arg}}})"
    else:
        latex_string = "0"

    print("\nThe remaining term of x is:")
    print(latex_string)


if __name__ == '__main__':
    solve_and_format()
