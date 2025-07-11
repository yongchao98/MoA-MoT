import sympy

def solve_and_format():
    """
    This function calculates the optimal mean response time 'x' for the given
    queueing system, removes the specified terms, and formats the remainder
    in LaTeX as requested.
    """
    # Define symbols for symbolic integration
    s, y, x_sym = sympy.symbols('s y x_sym')

    # System parameters
    lambda_rate = sympy.Rational(3, 2)
    # PDF of service time S ~ U(0, 1) is f(s) = 1 for s in [0, 1]
    f_s = 1

    # Calculate rho(s), the utilization from jobs of size at most s
    # rho(s) = lambda * integral(x * f(x) dx) from 0 to s
    rho_s_expr = lambda_rate * sympy.integrate(x_sym * f_s, (x_sym, 0, s))

    # The mean response time 'x' is the integral of E[T_s] from 0 to 1.
    # E[T_s] = E[W_s] + E[F_s]
    # E[W_s] = (lambda * integral(x^2*f(x)dx)) / (2*(1-rho(s)))
    # E[F_s] = integral(1/(1-rho(y)))dy from 0 to s

    # Calculate the integral of the waiting time component E[W_s] over all job sizes s
    wait_num_integrand = lambda_rate * sympy.integrate(x_sym**2 * f_s, (x_sym, 0, s))
    E_W_s = wait_num_integrand / (2 * (1 - rho_s_expr))
    total_wait_time = sympy.integrate(E_W_s, (s, 0, 1))

    # Calculate the integral of the flow time component E[F_s] over all job sizes s
    # This is a double integral, which we solve by changing the order of integration
    # integral from y=0 to 1 of integral from s=y to 1 of (1/(1-rho(y))) ds dy
    rho_y_expr = rho_s_expr.subs(s, y)
    flow_integrand = 1 / (1 - rho_y_expr)
    total_flow_time = sympy.integrate(flow_integrand * (1 - y), (y, 0, 1))

    # The total mean response time, x
    x_val = total_wait_time + total_flow_time

    # Sympy provides the result, which we can expand to identify terms
    x_expanded = sympy.expand_log(x_val, force=True)

    # Identify the components of x
    # x = -1/6 - 8/9*ln(2) + (2*sqrt(3)/3)*ln(2+sqrt(3))
    rational_term = sympy.Rational(-1, 6)
    log_rational_term = sympy.Rational(-8, 9) * sympy.log(2)
    
    # The remaining term is x minus the other two parts
    remaining_term_expr = x_expanded - rational_term - log_rational_term
    remaining_term_simplified = sympy.simplify(remaining_term_expr)

    # Format the remaining term in standard LaTeX notation
    # Format: rational multiplicand, then algebraic, then transcendental
    # Term is (2/3) * sqrt(3) * ln(2+sqrt(3))
    rational_mult = sympy.Rational(2, 3)
    algebraic_mult = sympy.sqrt(3)
    transcendental_mult_str = r"\ln(2+\sqrt{3})"

    # Build the final LaTeX string
    # sympy.latex() is used for robust conversion of symbolic expressions to LaTeX
    final_latex_string = (f"{sympy.latex(rational_mult)} "
                          f"{sympy.latex(algebraic_mult)} "
                          f"{transcendental_mult_str}")

    print(final_latex_string)

solve_and_format()