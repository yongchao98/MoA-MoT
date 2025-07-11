import sympy

def solve_sugra_parameters():
    """
    This script solves for the parameters alpha^2 and beta in the
    super-cosmological constant term of N=1, d=4 Supergravity.
    """

    # Define symbols for the parameters
    alpha, kappa, beta = sympy.symbols('alpha kappa beta')

    # --- Step 1: Set up the cancellation equations ---
    # The variation of the total Lagrangian is split into S-independent and S-linear parts.
    # For the action to be invariant, the non-total-derivative parts of the variation must cancel.

    # Equation from S-independent terms:
    # The variation contains a term (-1/(2*kappa) + alpha/4) * (e * epsilon_bar * gamma^mu * R_mu).
    # Since e*epsilon_bar*gamma^mu*R_mu is not a total derivative, its coefficient must vanish.
    # The other S-independent term, from the variation of the psi-bilinear, is already a total derivative.
    print("Step 1: Set up cancellation equations from SUSY variations.")
    print("--------------------------------------------------------")
    
    s_independent_eq = -1 / (2 * kappa) + alpha / 4
    print(f"S-independent cancellation equation: {s_independent_eq} = 0")

    # Equation from S-linear terms:
    # The coefficients of the term e*S*epsilon_bar*gamma^mu*psi_mu must sum to zero.
    # The terms come from delta(L_sugra), delta(alpha*e*S), and delta(alpha*e*kappa*beta*psi_bar*gamma*psi).
    # Their coefficients are -1/2, (3/4)*alpha*kappa, and alpha*kappa*beta respectively.
    s_linear_eq = -sympy.Rational(1, 2) + sympy.Rational(3, 4) * alpha * kappa + alpha * kappa * beta
    print(f"S-linear cancellation equation: {s_linear_eq} = 0")
    print("")

    # --- Step 2: Assume kappa = 1 (working in Planck units) ---
    print("Step 2: Assume kappa = 1 (Planck units).")
    print("-----------------------------------------")
    kappa_val = 1
    s_independent_eq = s_independent_eq.subs(kappa, kappa_val)
    s_linear_eq = s_linear_eq.subs(kappa, kappa_val)
    print(f"S-independent eq with kappa=1: {s_independent_eq} = 0")
    print(f"S-linear eq with kappa=1: {s_linear_eq} = 0")
    print("")

    # --- Step 3: Solve for alpha and alpha^2 ---
    print("Step 3: Solve for alpha and alpha^2.")
    print("------------------------------------")
    # Solve the S-independent equation for alpha
    alpha_sol = sympy.solve(s_independent_eq, alpha)
    if not alpha_sol:
        print("Could not solve for alpha.")
        return
    alpha_val = alpha_sol[0]
    alpha_sq_val = alpha_val**2
    print(f"From the S-independent equation, alpha = {alpha_val}")
    print(f"Therefore, alpha^2 = {alpha_sq_val}")
    print("")

    # --- Step 4: Solve for beta ---
    print("Step 4: Solve for beta.")
    print("-----------------------")
    # Substitute the value of alpha into the S-linear equation
    s_linear_eq_beta = s_linear_eq.subs(alpha, alpha_val)
    print(f"Substituting alpha = {alpha_val} into the S-linear equation:")
    print(f"{s_linear_eq_beta} = 0")
    
    # Solve for beta
    beta_sol = sympy.solve(s_linear_eq_beta, beta)
    if not beta_sol:
        print("Could not solve for beta.")
        return
    beta_val = beta_sol[0]
    print(f"Solving for beta, we get beta = {beta_val}")
    print("")

    # --- Final Results ---
    print("========================================")
    print("Final determined values of the parameters:")
    print(f"alpha^2 = {alpha_sq_val}")
    print(f"beta = {beta_val}")
    print("========================================")

if __name__ == '__main__':
    solve_sugra_parameters()