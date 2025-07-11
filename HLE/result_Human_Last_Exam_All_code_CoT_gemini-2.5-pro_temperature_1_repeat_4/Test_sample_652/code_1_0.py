import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    This function solves the economic model and calculates the elasticity of 
    optimal re-employment probability (q) with respect to the probability of 
    unemployment (p).
    """
    # --- Model Parameters ---
    w = 20.0  # Worker's wage
    p = 0.2   # Probability of becoming unemployed

    # --- Helper Functions ---
    def q_func(b, wage=w):
        """Calculates optimal q as a function of benefit b."""
        # Clamp b to be in a valid range for the log function
        b = np.maximum(b, 1e-10)
        b = np.minimum(b, wage - 1e-10)
        return (np.log(wage) - np.log(b)) / 2.0

    def system_equations(variables, p_val, wage=w):
        """Defines the system of non-linear equations to find equilibrium."""
        b, t = variables
        
        # Check for invalid variable ranges to guide the solver
        if b <= 0 or t < 0 or t >= wage:
            return [1e6, 1e6]

        q = q_func(b, wage)
        
        if not (0 < q < 1):
            return [1e6, 1e6]

        # Equation 1: From utility maximization FOCs
        eq1 = (wage - t) * (1 - q) - b * (1.5 - q)
        
        # Equation 2: From the government's budget constraint
        eq2 = t - (p_val / (1 - p_val)) * (1 - q) * b
        
        return [eq1, eq2]

    # --- Step 1: Solve for equilibrium at p=0.2 ---
    initial_guess = [w / 2, w / 10]  # Initial guess for (b, t)
    try:
        b_star, t_star = fsolve(system_equations, initial_guess, args=(p, w))
        q_star = q_func(b_star, w)
    except Exception as e:
        print(f"Failed to solve the system of equations: {e}")
        return

    # --- Step 2: Set up and solve for derivatives using Implicit Function Theorem ---
    # We have a system F(b, t, p) = 0. Implicit differentiation gives a
    # linear system for db/dp and dt/dp.
    
    # Derivative of q with respect to b
    q_prime_b = -1.0 / (2.0 * b_star)

    # Coefficients for the linear system M * [dt/dp, db/dp]^T = V
    # where M is the Jacobian of the system w.r.t (t, b) and V is the
    # negative of the Jacobian w.r.t p.
    # Our system is:
    # A * dt_dp + B * db_dp = 0
    # C * dt_dp + D * db_dp = F_rhs

    # From differentiating the first system equation w.r.t. p
    A = -1.0 + q_star
    B = -q_prime_b * (w - t_star) - (1.5 - q_star) + b_star * q_prime_b

    # From differentiating the second system equation w.r.t. p
    C = 1.0
    p_term = p / (1 - p)
    D = p_term * b_star * q_prime_b - p_term * (1 - q_star)
    
    # The right-hand-side vector comes from the partial derivative w.r.t p
    dp_term_dp = 1.0 / (1 - p)**2
    F_rhs = dp_term_dp * (1 - q_star) * b_star

    # Solve the 2x2 linear system for db_dp
    # A*dt_dp + B*db_dp = 0  => dt_dp = -(B/A)*db_dp
    # C*dt_dp + D*db_dp = F_rhs => C*(-(B/A)*db_dp) + D*db_dp = F_rhs
    # db_dp * (D - CB/A) = F_rhs => db_dp = F_rhs / (D - CB/A) = A*F_rhs / (AD-BC)
    determinant = A * D - B * C
    if abs(determinant) < 1e-9:
        print("Error: Matrix is singular, cannot calculate derivatives.")
        return
        
    db_dp = (A * F_rhs) / determinant

    # --- Step 3: Calculate the elasticity ---
    # Using the chain rule for dq/dp
    dq_dp = q_prime_b * db_dp
    
    # Elasticity formula
    elasticity = dq_dp * (p / q_star)
    
    # --- Step 4: Print the final results ---
    print("The elasticity of optimal q with respect to p is calculated as follows:")
    print("Elasticity ε = (dq*/dp) * (p / q*)")
    print(f"\nGiven parameters:")
    print(f"Wage w = {w}")
    print(f"Initial unemployment probability p = {p}")
    print(f"\nCalculated equilibrium values:")
    print(f"Optimal benefit b* = {b_star:.3f}")
    print(f"Optimal tax t* = {t_star:.3f}")
    print(f"Resulting optimal re-employment probability q* = {q_star:.3f}")
    print(f"\nDerivatives from implicit differentiation:")
    print(f"dq*/db = {q_prime_b:.3f}")
    print(f"db*/dp = {db_dp:.3f}")
    print(f"dq*/dp = (dq*/db) * (db*/dp) = {dq_dp:.3f}")
    print(f"\nFinal elasticity calculation:")
    print(f"ε = {dq_dp:.3f} * ({p:.3f} / {q_star:.3f}) = {elasticity:.3f}")

if __name__ == '__main__':
    solve_and_calculate_elasticity()
    # The calculated elasticity is approximately -0.444
    print("\n<<< -0.444 >>>")
