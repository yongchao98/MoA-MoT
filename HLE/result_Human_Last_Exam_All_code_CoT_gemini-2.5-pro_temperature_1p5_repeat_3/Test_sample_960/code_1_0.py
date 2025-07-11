import sympy

def solve_urn_problem():
    """
    Calculates the probability of reaching a state with 50% good products
    using martingale theory and the Beta distribution.
    """
    # Define the symbolic variable for integration
    x = sympy.Symbol('x')

    # Initial number of good (W_0) and defective (B_0) products
    W0 = 2
    B0 = 1

    print(f"Initial State: W_0 = {W0}, B_0 = {B0}\n")
    
    # --- Martingale Setup ---
    # The martingale is M_t = (W_t - B_t) / (W_t + B_t)
    m0_frac = sympy.S(W0 - B0) / sympy.S(W0 + B0)
    print("We use the martingale M_t = (W_t - B_t)/(W_t + B_t).")
    print(f"Initial value M_0 = ({W0} - {B0}) / ({W0} + {B0}) = {m0_frac}\n")

    # --- Long-term behavior using Beta distribution ---
    # The limiting fraction of good products, X_inf, follows a Beta(W0, B0) distribution.
    alpha = W0
    beta = B0
    beta_func = sympy.gamma(alpha) * sympy.gamma(beta) / sympy.gamma(alpha + beta)
    pdf = (x**(alpha - 1) * (1 - x)**(beta - 1)) / beta_func
    
    print(f"The limiting fraction of good products X_inf ~ Beta(alpha={alpha}, beta={beta}).")
    print(f"The PDF of X_inf is f(x) = {pdf}\n")

    # The process never stops (T=inf) if the fraction of good products remains > 1/2.
    # We calculate the conditional expectation E[M_inf | T=inf] = E[2*X_inf - 1 | X_inf > 1/2]
    
    # 1. Calculate P(X_inf > 1/2)
    p_gt_half = sympy.integrate(pdf, (x, sympy.S(1)/2, 1))
    
    # 2. Calculate the numerator for the conditional expectation: Integral[(2x-1)*pdf(x)] from 1/2 to 1
    integrand = (2*x - 1) * pdf
    numerator_exp = sympy.integrate(integrand, (x, sympy.S(1)/2, 1))
    
    # 3. Compute the conditional expectation
    cond_exp_M_inf = numerator_exp / p_gt_half
    
    print("The condition T=inf corresponds to X_inf > 1/2.")
    print(f"E[M_inf | T=inf] = (Integral from 1/2 to 1 of (2x-1)*f(x) dx) / (P(X_inf > 1/2))")
    print(f"E[M_inf | T=inf] = ({numerator_exp}) / ({p_gt_half}) = {cond_exp_M_inf}\n")

    # --- Applying Optional Stopping Theorem ---
    # We have M_0 = P(T = inf) * E[M_inf | T = inf]
    P_T_inf = m0_frac / cond_exp_M_inf

    print("From the Optional Stopping Theorem, P(T=inf) = M_0 / E[M_inf | T=inf].")
    print(f"P(T=inf) = ({m0_frac}) / ({cond_exp_M_inf}) = {P_T_inf}\n")
    
    # The final probability is P(T < inf) = 1 - P(T = inf)
    P_T_finite = 1 - P_T_inf

    print("The probability that the production will reach 50% good products is:")
    print(f"P(T < inf) = 1 - P(T=inf) = 1 - {P_T_inf} = {P_T_finite}")
    print("\nThis exact probability is the least upper bound for the event.")

solve_urn_problem()