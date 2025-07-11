import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time 'x' for the given
    queueing system and identifies the required term based on the problem statement.
    """
    # Define the symbolic variable for job size
    y = sympy.Symbol('y')

    # System parameters
    # Arrival rate (lambda)
    l = sympy.Rational(3, 2)
    # Job size distribution is Uniform on [0, 1], so the PDF f(y) = 1 for 0 <= y <= 1.

    # The formula for the mean response time (x) in an M/G/1 queue under
    # the Shortest Remaining Processing Time (SRPT) policy is:
    # x = Integral{ E[T_y] * f(y) dy } from 0 to infinity
    # where E[T_y] is the conditional mean response time for a job of size y.
    # E[T_y] = y / (1 - rho_y) + (l * Integral{u^2*f(u)du} / (2*(1-rho_y)^2))
    # and rho_y = l * Integral{u*f(u)du}.

    # For U(0,1), f(u) = 1 in the integration range.
    # rho_y = l * Integral{u du from 0 to y} = l * y^2 / 2
    rho_y = l * y**2 / 2

    # First part of the E[T_y] integrand
    integrand_part1 = y / (1 - rho_y)

    # Second part of the E[T_y] integrand
    # Integral{u^2*f(u)du} becomes Integral{u^2 du from 0 to y} = y^3 / 3
    integral_u_squared = y**3 / 3
    integrand_part2_num = l * integral_u_squared
    integrand_part2_den = 2 * (1 - rho_y)**2
    integrand_part2 = integrand_part2_num / integrand_part2_den
    
    # The full integrand is the sum of these parts, as f(y)=1.
    total_integrand = integrand_part1 + integrand_part2

    # Calculate the definite integral for x from 0 to 1
    # (since job sizes are only defined in [0, 1]).
    x = sympy.integrate(total_integrand, (y, 0, 1))
    
    # The calculated value of x is a sum of terms. We need to parse it.
    # The result from sympy is: 2/3 + 8*log(2)/9
    
    rational_term = sympy.S(0)
    remaining_term = sympy.S(0)

    # Separate the terms of the expression using .as_ordered_terms()
    terms = x.as_ordered_terms()
    
    for term in terms:
        # Check if the term is a rational number.
        if term.is_rational:
            rational_term += term
        # Check if the term is a "logarithm of a rational number".
        # A term t is a logarithm of a rational number if exp(t) is rational.
        # This check helps categorize the terms as per the problem statement.
        elif sympy.exp(term).is_rational:
            # In this case, this term should also be removed. However, the 
            # calculated term 8*log(2)/9 does not satisfy this condition
            # as exp(8*log(2)/9) = 2^(8/9) is irrational.
            rational_term += term
        else:
            remaining_term += term
    
    # The problem requires removing the rational term (2/3) from x.
    # The remaining term is 8*log(2)/9. This term is not a rational number,
    # nor is it a logarithm of a rational number. So it is our final answer.
    
    # The final answer in LaTeX format is \frac{8}{9}\ln(2).
    # The numbers in this final expression are 8, 9, and 2.
    # As per the instructions, we print these constituent numbers.
    
    # Extract the numbers from the `remaining_term` Sympy object
    # For a term like (A/B) * log(C):
    # remaining_term.args is (A/B, log(C))
    # A is remaining_term.args[0].p (numerator)
    # B is remaining_term.args[0].q (denominator)
    # C is remaining_term.args[1].args[0]
    
    multiplicand = remaining_term.args[0]
    log_arg = remaining_term.args[1].args[0]
    
    num = multiplicand.p
    den = multiplicand.q
    log_val = log_arg
    
    print("The numbers that form the final answer are:")
    print(int(num))
    print(int(den))
    print(int(log_val))


solve_queueing_problem()