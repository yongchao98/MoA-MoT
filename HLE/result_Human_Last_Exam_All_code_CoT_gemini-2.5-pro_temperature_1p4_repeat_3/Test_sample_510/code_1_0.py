import sympy as sp

def solve_queueing_problem():
    """
    Calculates the optimal mean response time 'x' for the given M/G/1 SRPT queue
    and isolates the required term based on the problem's criteria.
    """
    # Define symbols for symbolic mathematics
    s, t = sp.symbols('s t', real=True, positive=True)

    # --- Step 1: Define system parameters ---
    # Arrival rate
    lambda_val = sp.Rational(3, 2)
    # Service time distribution is U(0, 1), so its PDF is 1 on the interval [0, 1]
    f_s = 1
    # The integration limit for job sizes
    service_time_upper_bound = 1

    # --- Step 2: Calculate Mean Service Time E[S] ---
    E_S = sp.integrate(s * f_s, (s, 0, service_time_upper_bound))

    # --- Step 3: Calculate terms for Mean Waiting Time E[W] ---
    # ρ_s: traffic intensity from jobs of size at most s
    rho_s = lambda_val * sp.integrate(t * f_s, (t, 0, s))

    # Numerator term for E[V_s]: λ * ∫[0 to s] t^2 f(t) dt
    numerator_term_V_s = lambda_val * sp.integrate(t**2 * f_s, (t, 0, s))
    
    # E[V_s]: mean waiting time for a job of size s
    E_V_s = numerator_term_V_s / (2 * (1 - rho_s))

    # Integrand for E[W]
    integrand_W = E_V_s * f_s

    # --- Step 4: Calculate Mean Waiting Time E[W] ---
    E_W = sp.integrate(integrand_W, (s, 0, service_time_upper_bound))

    # --- Step 5: Calculate optimal mean response time x = E[T] ---
    x = E_S + E_W

    print(f"The optimal mean response time is x = E[S] + E[W].")
    print(f"Calculated E[S] = {E_S}")
    print(f"Calculated E[W] = {E_W}")
    print(f"Therefore, x = {x}\n")

    # --- Step 6: Isolate the required term from x ---
    # The value of x is a sum of terms. We need to remove additive rational terms
    # and additive terms which are logarithms of rational numbers.
    terms = x.as_ordered_terms()
    
    rational_part = 0
    remaining_term = 0
    
    for term in terms:
        # Identify and separate the purely rational part
        if term.is_rational:
            rational_part += term
        else:
            # According to the problem, we must check if this term is a logarithm
            # of a rational number. A term 'y' is a log of a rational number if 
            # exp(y) is rational.
            # For the term 4/9*ln(2), exp(4/9*ln(2)) = 2^(4/9), which is not rational.
            # Thus, this term is not removed.
            remaining_term += term
    
    print("The expression for x can be decomposed into:")
    print(f"  - An additive rational term: {rational_part}")
    print(f"  - A remaining term: {remaining_term}\n")
    print("After removing the additive rational term, we are left with the final answer.")
    
    # --- Step 7: Format the final answer ---
    # The problem asks for the equation and the final term in LaTeX.
    # Extract the numbers from the remaining term to construct the "equation".
    if remaining_term:
        coeff, trans_part = remaining_term.as_coeff_Mul()
        num = coeff.p
        den = coeff.q
        log_arg = trans_part.args[0]
    
        # The final equation showing the numbers involved
        print("Final Equation for the Remaining Term:")
        print(f"({num}/{den}) * ln({log_arg})")

        # The final answer in LaTeX format
        latex_answer = sp.latex(remaining_term)
        print("\nIn standard LaTeX notation, this is:")
        print(latex_answer)


if __name__ == '__main__':
    solve_queueing_problem()
