import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time x and processes it
    according to the problem's instructions.
    """
    # 1. Define symbols and constants
    s, y = sympy.symbols('s y')
    lambda_val = sympy.Rational(3, 2)

    # 2. Define rho(s), the traffic intensity from jobs of size <= s
    # f(y) = 1 for y in [0,1]
    rho_s = lambda_val * sympy.integrate(y, (y, 0, s))

    # 3. Define the two parts of the integrand for E[T(s)]
    # First part corresponds to s / (1 - rho(s))
    integrand1 = s / (1 - rho_s)
    
    # Second part corresponds to the waiting time component
    numerator_term = lambda_val * sympy.integrate(y**2, (y, 0, s))
    integrand2 = numerator_term / (2 * (1 - rho_s)**2)

    # 4. Calculate the two integrals that form x
    I1 = sympy.integrate(integrand1, (s, 0, 1))
    I2 = sympy.integrate(integrand2, (s, 0, 1))
    
    # The optimal mean response time is the sum
    x = I1 + I2

    # 5. Output the final equation and its components as requested
    print("The optimal mean response time, x, is the sum of two integrals, I1 and I2.")
    print(f"I1 = Integral from 0 to 1 of ({sympy.simplify(integrand1)}) ds = {I1}")
    print(f"I2 = Integral from 0 to 1 of ({sympy.simplify(integrand2)}) ds = {I2}")
    
    print("\nThe final equation for x is:")
    # Using sstr for a compact string representation of the SymPy expressions
    i1_str = sympy.sstr(I1, full_prec=False)
    i2_str = sympy.sstr(I2, full_prec=False)
    x_str = sympy.sstr(x, full_prec=False)
    print(f"x = I1 + I2 = ({i1_str}) + ({i2_str})")
    print(f"x = {x_str}")

    # 6. Analyze x and remove terms as instructed
    print("\nRemoving terms from x:")
    # x is composed of additive terms. We iterate through them.
    terms = sympy.Add.make_args(x)
    remaining_term = sympy.S(0)

    for term in terms:
        # Check if the term is rational
        if term.is_rational:
            print(f"- Removing rational term: {term}")
            continue
            
        # Check if the term T is a logarithm of a rational number, i.e., T = ln(q) for q rational.
        # This is equivalent to checking if exp(T) is rational.
        is_log_of_rational = False
        try:
            exp_term = sympy.exp(term)
            if exp_term.is_rational:
                is_log_of_rational = True
        except TypeError:
            # sympy may raise an error for certain complex expressions
            pass

        if is_log_of_rational:
            print(f"- Removing term which is a logarithm of a rational number: {term}")
            continue
        else:
            print(f"- Keeping term: {term}, because it is not a rational number nor a logarithm of a rational number.")
            exp_of_term = sympy.exp(term)
            print(f"  (exp({term}) = {exp_of_term}, which is not rational)")

        # If not removed, it's part of the answer
        remaining_term += term

    # 7. Format the final answer in LaTeX
    final_answer_latex = sympy.latex(remaining_term, ln_notation=True)
    
    print(f"\nThe remaining term of x is: {remaining_term}")
    print("\nFinal Answer in LaTeX format:")
    print(f"<<<{final_answer_latex}>>>")

if __name__ == '__main__':
    solve_queueing_problem()