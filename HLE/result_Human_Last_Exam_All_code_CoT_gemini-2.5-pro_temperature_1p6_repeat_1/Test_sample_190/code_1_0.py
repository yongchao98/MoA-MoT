import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.

    The method uses a diffusion approximation for the Markov chain for large states k.
    1.  Calculates the drift (mu_k) and the second moment of jumps (sigma_k_sq).
    2.  Sets up the corresponding drift (b(x)) and diffusion (D(x)) coefficients
        for the approximating diffusion process.
    3.  Applies the transience condition for a diffusion on [a, infinity) with a
        reflecting boundary at a. Transience occurs if the process escapes to
        infinity, which depends on the convergence of an integral involving
        b(x) and D(x).
    4.  This convergence condition yields an inequality for c. The infimum of the
        set of solutions for c is the answer.
    """
    # Define symbols
    k, c = sympy.symbols('k c')

    # Transition probabilities for a jump of size j from state k
    # P_k,k+j
    transitions = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
         1: sympy.Rational(1, 4) + c/k,
         2: sympy.Rational(1, 4)
    }

    # Step 1: Calculate the drift (mu_k) and the second moment (sigma_k_sq)
    mu_k = sympy.simplify(sum(jump * prob for jump, prob in transitions.items()))
    sigma_k_sq = sympy.simplify(sum(jump**2 * prob for jump, prob in transitions.items()))

    # Step 2: Define coefficients for the diffusion approximation (for variable x)
    x = sympy.Symbol('x')
    b_x = mu_k.subs(k, x)
    D_x = (sigma_k_sq / 2)

    # Step 3: Analyze transience.
    # The scale density s'(x) is proportional to exp(-integral(b(y)/D(y) dy)).
    # The process is transient (escapes to oo) if integral of s'(x) from a to oo converges.
    # This integral converges if the exponent of x in s'(x) is < -1.
    y = sympy.Symbol('y')
    ratio_b_D = sympy.simplify(b_x.subs(x, y) / D_x)
    
    # The exponent in the scale density is the integral of the ratio
    exponent_integral = sympy.integrate(ratio_b_D, y)
    
    # The scale density is proportional to x raised to some power. Find that power.
    # The scale density is exp(-exponent_integral)
    scale_density_power = -sympy.Poly(exponent_integral, sympy.log(y)).coeffs()[0]

    # Step 4: Apply the convergence criterion for the integral of the scale density.
    # The integral of x^p from a to oo converges if p < -1.
    # Here p is the scale_density_power.
    inequality = sympy.solve(scale_density_power < -1, c)

    # Step 5: Determine the infimum from the resulting inequality
    # The inequality will be of the form c > value. The infimum is that value.
    if inequality.rel_op == '>':
        infimum_c = inequality.rhs
    elif inequality.rel_op == '<':
        infimum_c = inequality.lhs
    else: # Should not happen for this problem
        infimum_c = None

    print("Analysis of the Markov Chain Transience:")
    print("------------------------------------------")
    print("1. Asymptotic properties for large state k:")
    print(f"   - Expected jump (drift) mu_k = {mu_k}")
    print(f"   - Expected squared jump sigma_k^2 = {sigma_k_sq}")
    print("\n2. Diffusion Approximation:")
    print(f"   - Drift coefficient b(x) = {b_x}")
    print(f"   - Diffusion coefficient D(x) = sigma_x^2 / 2 = {D_x}")
    print("\n3. Transience Condition:")
    print("   The chain is transient if the approximating diffusion process escapes to infinity.")
    print("   This occurs if the integral of the scale density s'(x) from a to infinity converges.")
    print(f"   The scale density s'(x) is proportional to x^p where p = {scale_density_power}.")
    print("\n4. Deriving the inequality for c:")
    print("   The integral converges if the power p is less than -1.")
    print(f"   Solving the inequality {scale_density_power} < -1 for c:")
    print(f"   {scale_density_power.lhs}*c < -1")
    print(f"   -8*c < -5")
    print(f"   8*c > 5")
    print(f"   c > 5/8")
    print(f"   So, the chain is transient for {inequality}.")

    print("\n5. Conclusion:")
    print("   The set of values for c where the chain is transient is (5/8, infinity).")
    print("   The infimum (greatest lower bound) of this set is 5/8.")

    final_c_num = infimum_c.p
    final_c_den = infimum_c.q

    print("\nFinal Answer:")
    print("The final equation for the infimum c is:")
    print(f"c = {final_c_num} / {final_c_den}")


solve_markov_chain_transience()