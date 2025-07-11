import sympy

def solve_queueing_problem():
    """
    Calculates the optimal mean response time x for the given single-server queueing system,
    then isolates and formats the required remaining term according to the problem description.
    """
    # 1. Define symbols and parameters
    s, y = sympy.symbols('s y', real=True, positive=True)
    lam = sympy.Rational(3, 2)  # Arrival rate lambda
    # Service time PDF f(s) is 1 for s in (0, 1), and 0 otherwise.
    # The integration bounds (0, 1) will account for this.

    # 2. Define rho(s), the traffic intensity from jobs of size <= s
    # rho(s) = lambda * integral from 0 to s of y*f(y) dy
    # Since f(y)=1, this is lambda * integral of y
    rho_s = lam * sympy.integrate(y, (y, 0, s))

    # 3. Define the integrands for the two parts of the mean response time formula
    # Part 1: I1
    integrand1 = s / (1 - rho_s)

    # Part 2: I2
    # Inner integral for I2: integral from 0 to s of y^2*f(y) dy
    inner_integral_I2 = sympy.integrate(y**2, (y, 0, s))
    integrand2 = (lam * inner_integral_I2) / (2 * (1 - rho_s)**2)

    # 4. Calculate the definite integrals from 0 to 1
    I1 = sympy.integrate(integrand1, (s, 0, 1))
    I2 = sympy.integrate(integrand2, (s, 0, 1))

    # 5. The optimal mean response time x is the sum of the two integrals
    x = I1 + I2

    # Print the intermediate steps as requested
    print("The optimal mean response time, x, is the sum of two integrals, I1 and I2.")
    print(f"The first integral I1 evaluates to: {sympy.latex(I1)}")
    print(f"The second integral I2 evaluates to: {sympy.latex(I2)}")
    print(f"The full expression for x is: x = I1 + I2 = {sympy.latex(x)}")
    print("-" * 20)

    # 6. Process x to find the remaining term
    # The expression for x is a sum of terms. We iterate through them.
    remaining_terms = []
    for term in x.as_ordered_terms():
        # Check if the term is a rational number
        if term.is_rational:
            print(f"Removing additive rational term: {sympy.latex(term)}")
            continue

        # Check if the term is a logarithm of a rational number.
        # A term 't' is a log of a rational number if t = log(q) where q is rational.
        # This is equivalent to checking if exp(t) is rational.
        is_log_of_rational = False
        try:
            if sympy.exp(term).is_rational:
                is_log_of_rational = True
        except TypeError: # sympy might not be able to determine rationality for complex expressions
            pass
        
        if is_log_of_rational:
            print(f"Removing term which is a logarithm of a rational number: {sympy.latex(term)}")
            continue

        # If the term is not removed, add it to our list of remaining terms
        remaining_terms.append(term)

    # 7. Format the final answer
    final_answer_expr = sum(remaining_terms)
    # Convert to a prettier LaTeX string
    final_answer_latex = sympy.latex(final_answer_expr, ln_notation=True)
    # The default latex for a/b * ln(c) is \frac{a \ln{c}}{b}, which is good.
    # Let's make it match the requested format exactly.
    coeff, rest = final_answer_expr.as_coeff_mul()
    if rest and rest[0].func == sympy.log:
        log_arg = rest[0].args[0]
        final_answer_latex = f"{sympy.latex(coeff)} \\ln {sympy.latex(log_arg)}"

    print(f"The remaining term is: {final_answer_latex}")
    print(f"<<<{final_answer_latex}>>>")

if __name__ == '__main__':
    solve_queueing_problem()