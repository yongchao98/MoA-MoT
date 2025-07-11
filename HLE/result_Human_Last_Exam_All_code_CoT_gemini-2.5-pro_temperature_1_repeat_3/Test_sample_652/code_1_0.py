import sympy as sp

def solve_elasticity():
    """
    This function solves the economic problem to find the elasticity of optimal q with respect to p.
    """
    # Step 1: Define symbols and fundamental relationships
    # b: unemployment benefit
    # p: probability of becoming unemployed
    b, p = sp.symbols('b p', real=True, positive=True)
    w = 20  # Wage
    p_val = 0.2 # Given probability

    # The worker chooses q to maximize q*ln(w) + (1-q)*ln(b) - q^2
    # FOC: ln(w) - ln(b) - 2q = 0
    # This gives optimal q as a function of b
    q_star_expr = (sp.ln(w) - sp.ln(b)) / 2

    # Step 2: Define the government's budget constraint
    # tax revenue (t) = expected benefit payments
    t_expr = p * b * (1 - q_star_expr)

    # Step 3: Define the government's objective function (worker's overall expected utility)
    # We can simplify the utility in the unemployed state: U_u = (q*)^2 + ln(b)
    unemployed_utility_expr = q_star_expr**2 + sp.ln(b)
    # Overall expected utility
    EU_expr = (1 - p) * sp.ln(w - t_expr) + p * unemployed_utility_expr

    # Step 4: Find the optimal benefit b* by solving the FOC
    # Differentiate the expected utility w.r.t. b to get the FOC for the government's problem
    FOC_expr = sp.diff(EU_expr, b)

    # Substitute p = 0.2 and solve for b numerically.
    # The search space for b is where 0 < q < 1, which implies 20/e^2 < b < 20.
    # A starting guess of 10 is appropriate.
    FOC_at_p = FOC_expr.subs(p, p_val)
    try:
        b_star_val = sp.nsolve(FOC_at_p, b, 10)
    except ValueError:
        print("Could not find a numerical solution for b*.")
        return

    # Step 5: Calculate the elasticity
    # Elasticity E = (dq*/dp) * (p / q*)

    # First, calculate q* at the optimal b
    q_star_val = q_star_expr.subs(b, b_star_val)

    # To find dq*/dp, we use the chain rule: dq*/dp = (dq*/db) * (db*/dp)
    # We find db*/dp using the implicit function theorem on the FOC:
    # Let G(b, p) = FOC_expr. Then db/dp = - (dG/dp) / (dG/db)
    
    dG_dp_expr = sp.diff(FOC_expr, p)
    dG_db_expr = sp.diff(FOC_expr, b) # This is the second-order condition

    # Evaluate these derivatives at the optimal point (b*, p_val)
    subs_dict = {p: p_val, b: b_star_val}
    dG_dp_val = dG_dp_expr.subs(subs_dict)
    dG_db_val = dG_db_expr.subs(subs_dict)

    # db*/dp
    db_dp_val = -dG_dp_val / dG_db_val

    # dq*/db
    dq_db_expr = sp.diff(q_star_expr, b)
    dq_db_val = dq_db_expr.subs(b, b_star_val)

    # dq*/dp
    dq_dp_val = dq_db_val * db_dp_val
    
    # Finally, the elasticity
    elasticity = dq_dp_val * (p_val / q_star_val)

    # Print the final result rounded to three decimal places
    print(f"{elasticity:.3f}")

if __name__ == "__main__":
    solve_elasticity()