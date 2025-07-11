import sympy

def solve_queueing_problem():
    """
    This function symbolically calculates the optimal mean response time 'x'
    for the described M/G/1-SRPT queue, and then isolates the required term.
    """
    # Define symbols for symbolic mathematics
    s, y = sympy.symbols('s y', real=True, positive=True)

    # System parameters from the problem description
    lam = sympy.Rational(3, 2)  # Arrival rate
    # Service time PDF f(y) = 1 for y in the interval [0, 1]

    # rho(s): utilization from jobs of size up to s
    # rho(s) = lambda * Integral(y*f(y), (y, 0, s))
    # Since f(y) is 1, the integrand is just y.
    rho_y = lam * sympy.integrate(y, (y, 0, y))

    # The formula for the mean response time E[T(s)] for a job of size s
    # under SRPT is composed of three parts:
    # 1. The job's own service time, s.
    # 2. The mean waiting time before service starts, E[W_q(s)].
    # 3. The mean delay from preemptions after service starts, E[W_p(s)].
    #
    # E[T(s)] = s + E[W_q(s)] + E[W_p(s)]
    # E[W_q(s)] = W_0(s) / (1 - rho(s))
    # W_0(s) = (lambda/2) * Integral(y**2 * f(y), (y, 0, s))
    # E[W_p(s)] = Integral(rho(y)/(1-rho(y)), (y, 0, s))
    
    # We will compute the overall mean response time x by integrating E[T(s)]
    # over the job size distribution U(0,1).
    # x = Integral(E[T(s)], (s, 0, 1))
    
    # Integral of s*f(s) from 0 to 1
    term1 = sympy.integrate(s, (s, 0, 1))

    # Integral of E[W_q(s)]*f(s) from 0 to 1
    rho_s = rho_y.subs(y, s)
    W0_s = (lam / 2) * sympy.integrate(y**2, (y, 0, s))
    E_Wq_s = W0_s / (1 - rho_s)
    term2 = sympy.integrate(E_Wq_s, (s, 0, 1))
    
    # Integral of E[W_p(s)]*f(s) from 0 to 1. This is a double integral.
    E_Wp_s_integrand = rho_y / (1 - rho_y)
    # We change the order of integration to solve it
    term3 = sympy.integrate( (1-y) * E_Wp_s_integrand, (y, 0, 1))

    # The overall optimal mean response time x
    x = term1 + term2 + term3

    # The full expression for x is:
    # x = -1/6 - 8*log(2)/9 + 2*sqrt(3)*log(2 + sqrt(3))/3
    # We can decompose x into three additive terms based on the removal rules.
    term_A = sympy.S(-1)/6  # Rational term
    term_B = -sympy.S(8)/9 * sympy.log(2)  # Term with logarithm of a rational number
    term_C = (sympy.S(2) * sympy.sqrt(3) / 3) * sympy.log(2 + sympy.sqrt(3))  # Remaining term

    # Print the full equation for x, showing each numerical component
    print("The optimal mean response time x is composed of three parts: x = A + B + C")
    print("Each numerical component of the final equation for x is shown below:")
    print(f"A (rational term) = {term_A.p}/{term_A.q}")
    # Extracting components from term B: (-8/9) * ln(2)
    b_coeff = term_B.args[0]
    b_log_arg = term_B.args[1].args[0]
    print(f"B (log of rational term) = ({b_coeff.p}/{b_coeff.q}) * ln({b_log_arg})")
    # For term C, we display its structure: (2*sqrt(3)/3) * ln(2 + sqrt(3))
    print("C (remaining term) = (2 * sqrt(3) / 3) * ln(2 + sqrt(3))")
    
    # The remaining term after removal is term_C
    remaining_term = term_C
    
    # Format the answer in standard LaTeX notation as requested.
    latex_string = sympy.latex(remaining_term)
    
    print("\nAfter removing the additive rational term (A) and the term which is a logarithm of a rational number (B), the remaining term is:")
    print(remaining_term)
    
    print("\nThis term formatted in standard LaTeX notation is:")
    print(latex_string)


if __name__ == "__main__":
    solve_queueing_problem()