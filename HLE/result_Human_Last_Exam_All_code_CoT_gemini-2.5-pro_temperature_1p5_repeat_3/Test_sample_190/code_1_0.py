import sympy

def solve_transience_problem():
    """
    This script calculates the infimum of c for which the given Markov chain is transient.
    The steps are explained through print statements.
    """

    # Step 1: Define symbolic variables
    c = sympy.Symbol('c')
    k = sympy.Symbol('k')

    # Step 2: Calculate the expected drift (mu_k) for large k.
    # The jump values are: -2, +2, -1, +1
    # The probabilities are: 1/4, 1/4, 1/4 - c/k, 1/4 + c/k
    
    jump_values = [-2, 2, -1, 1]
    probabilities = [sympy.Rational(1, 4), sympy.Rational(1, 4), 
                     sympy.Rational(1, 4) - c/k, sympy.Rational(1, 4) + c/k]

    mu_k = sum(jump * prob for jump, prob in zip(jump_values, probabilities))
    mu_k = sympy.simplify(mu_k)

    print("Step 1: Calculate the expected drift mu(k) for large k.")
    print("mu(k) = E[X_{n+1} - X_n | X_n = k]")
    print(f"mu(k) = ({jump_values[0]})*({probabilities[0]}) + ({jump_values[1]})*({probabilities[1]}) + ({jump_values[2]})*({probabilities[2]}) + ({jump_values[3]})*({probabilities[3]})")
    print(f"mu(k) = {mu_k}")
    print("-" * 30)

    # Step 3: Calculate the expected squared jump distance (m2_k)
    m2_k = sum(jump**2 * prob for jump, prob in zip(jump_values, probabilities))
    m2_k = sympy.simplify(m2_k)
    
    print("Step 2: Calculate the expected squared jump distance m2(k) for large k.")
    print("m2(k) = E[(X_{n+1} - X_n)^2 | X_n = k]")
    print(f"m2(k) = ({jump_values[0]})**2*({probabilities[0]}) + ({jump_values[1]})**2*({probabilities[1]}) + ({jump_values[2]})**2*({probabilities[2]}) + ({jump_values[3]})**2*({probabilities[3]})")
    print(f"m2(k) = {m2_k}")
    print("-" * 30)

    # Step 4: Apply the transience criterion
    # The chain is transient if the sum of exp(-Integral(2*mu(x)/m2(x) dx)) over k converges.
    # The term inside the integral is:
    exponent_integrand = 2 * mu_k / m2_k
    
    print("Step 3: Apply the Lamperti/Pakes criterion for transience.")
    print("The chain is transient if Sum_k exp(-Integral(2*mu(x)/m2(x) dx)) converges.")
    print("Let's first compute the term to be integrated:")
    print(f"2 * mu(k) / m2(k) = 2 * ({mu_k}) / ({m2_k}) = {exponent_integrand}")
    
    # Integrate with respect to k
    exponent_integral = sympy.integrate(exponent_integrand, k)
    
    print("The integral is:")
    print(f"Integral({exponent_integrand}, k) = {exponent_integral}")
    
    # The term in the sum behaves like k^(-p)
    # where p is the coefficient of log(k)
    p = exponent_integral.coeff(sympy.log(k))
    print(f"The term to be summed is approximately k^(-p), where p = {p}")
    print("-" * 30)
    
    # Step 5: Analyze the convergence of the resulting p-series.
    print("Step 4: Determine the condition on c for the series to converge.")
    print(f"The series Sum_k k^(-p) converges if p > 1.")
    # The numbers in the final equation are p's numerator and denominator, and 1.
    p_num, p_den = p.as_numer_denom()
    print(f"The condition for transience is: {p_num}*c / {p_den} > 1")
    print(f"Solving for c: {p_num}*c > {p_den}")
    c_inequality_rhs_num = p_den
    c_inequality_rhs_den = p_num
    print(f"==> c > {c_inequality_rhs_num}/{c_inequality_rhs_den}")
    print("-" * 30)

    # Step 6: Find the infimum.
    infimum = sympy.Rational(c_inequality_rhs_num, c_inequality_rhs_den)
    print("Step 5: Find the infimum.")
    print(f"The set of c for which the chain is transient is ( {infimum}, oo ).")
    print(f"The infimum of this set is {infimum}.")
    print(f"The final numerical answer is: {float(infimum)}")

solve_transience_problem()
<<<0.625>>>