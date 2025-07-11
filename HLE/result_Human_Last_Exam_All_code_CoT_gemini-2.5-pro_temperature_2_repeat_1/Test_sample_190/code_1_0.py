import sympy

def find_transience_infimum():
    """
    This script finds the infimum of the set of values 'c' for which a given
    Markov chain is transient.

    The method involves:
    1. Calculating the first and second moments of the jumps from a state k.
    2. Setting up a differential equation for a Lyapunov function f(k) based on
       these moments.
    3. Solving the ODE to find the form of f(k).
    4. Applying Lamperti's criterion for transience/recurrence based on whether
       f(k) is bounded or unbounded as k -> infinity.
    5. Solving for the critical value of c that separates these two behaviors.
    """
    # Define symbols for state k and parameter c
    k, c = sympy.symbols('k c', real=True, positive=True)

    # The transition probabilities from state k are given for jumps of -2, -1, 1, 2.
    jumps_and_probs = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
        1:  sympy.Rational(1, 4) + c/k,
        2:  sympy.Rational(1, 4)
    }

    # Step 1: Calculate the first moment (drift), mu_k.
    mu_k = sympy.simplify(sum(jump * prob for jump, prob in jumps_and_probs.items()))
    print("Step 1: Calculate the drift mu_k = E[X_{n+1} - k | X_n=k]")
    print("mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print(f"mu_k = {mu_k}\n")

    # Step 2: Calculate the second moment, M_2k.
    M_2k = sympy.simplify(sum(jump**2 * prob for jump, prob in jumps_and_probs.items()))
    print("Step 2: Calculate the second moment M_2k = E[(X_{n+1} - k)^2 | X_n=k]")
    print("M_2k = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)")
    print(f"M_2k = {M_2k}\n")

    # Step 3: Set up the differential equation for the test function f(k).
    # The condition E[f(X_{n+1}) - f(k)] ≈ 0 translates to f'(k)*mu_k + (1/2)*f''(k)*M_2k = 0.
    # Let g(k) = f'(k), so g'(k) = f''(k). The equation becomes g(k)*mu_k + (1/2)*g'(k)*M_2k = 0.
    g = sympy.Function('g')
    ode = (M_2k / 2) * g(k).diff(k) + mu_k * g(k)
    print("Step 3: Set up the ODE for g(k) = f'(k): (M_2k/2)*g'(k) + mu_k*g(k) = 0")
    print(f"Substituting the moments, the equation is: {ode} = 0\n")

    # Step 4: Solve the ODE to find the form of g(k) = f'(k).
    # The solution determines the integrand for f(k).
    # The equation is: g'(k)/g(k) = -2*mu_k / M_2k
    exponent = sympy.simplify(-2 * mu_k / M_2k * k) # We solve d(ln(g))/dk = -2*mu_k/M_2k
    print(f"Step 4: Solve the ODE. The solution for g(k) is of the form C*k^p.")
    print(f"The exponent p is given by k * (-2 * mu_k / M_2k) = k * (-2 * ({mu_k}) / ({M_2k})) = {exponent}")
    print(f"So, f'(k) is proportional to k^({exponent}).\n")
    
    # Step 5: Apply Lamperti's criterion and find the critical value for c.
    # The function f(k) is the integral of f'(k). f(k) ~ integral(k^p dk) ~ k^(p+1).
    # For transience, f(k) must be bounded as k -> infinity, which means the exponent p+1 < 0.
    # For recurrence, f(k) must be unbounded, which means p+1 >= 0.
    # The boundary case is when the exponent of k in f(k) is 0, i.e., p + 1 = 0.
    final_exponent = exponent + 1
    boundary_eq = sympy.Eq(final_exponent, 0)
    
    print("Step 5: Find the boundary for c using Lamperti's criterion.")
    print("The chain is transient if f(k) is bounded, which means the exponent in f(k) must be < 0.")
    print(f"The exponent of k in f(k) is {exponent} + 1 = {final_exponent}.")
    print(f"The boundary condition for transience vs. recurrence is when this exponent equals 0:")
    print(f"{final_exponent} = 0")
    
    # Solve for c at the boundary
    critical_c_sol = sympy.solve(boundary_eq, c)
    critical_c = critical_c_sol[0]
    
    print(f"\nSolving the equation for c yields c = {critical_c}.")
    print(f"The Markov chain is transient if c > {critical_c} and recurrent if c <= {critical_c}.")
    print(f"Thus, the set of c for which the chain is transient is ({critical_c}, oo).")
    print(f"The infimum of this set is {critical_c}.")

if __name__ == '__main__':
    find_transience_infimum()
    print("\n<<<5/8>>>")