import numpy as np
from scipy.optimize import minimize

def solve_race_growth():
    """
    Calculates the difference between the optimal (W*) and actual (W) 
    growth rates for a horse race betting problem with mistaken probabilities.
    """
    
    # Define payout odds (b to 1)
    # b = [b_A, b_B, b_C]
    b = np.array([4.0, 3.0, 3.0])

    # Define true probabilities
    p_true = np.array([0.5, 0.25, 0.25])
    
    # Define mistaken (believed) probabilities
    p_believed = np.array([0.25, 0.5, 0.25])

    # Define the objective function to minimize (negative of expected log growth)
    def neg_log_growth(f, p, b):
        """
        Calculates the negative of the expected log growth rate.
        f: fractions of bankroll bet [f_A, f_B, f_C]
        p: probabilities of winning [p_A, p_B, p_C]
        b: net odds [b_A, b_B, b_C]
        """
        # Ensure fractions are handled as a numpy array
        f = np.asarray(f)
        
        # Calculate wealth growth factors for each outcome
        g_A = 1 - f[1] - f[2] + b[0] * f[0]
        g_B = 1 - f[0] - f[2] + b[1] * f[1]
        g_C = 1 - f[0] - f[1] + b[2] * f[2]
        
        # If any growth factor is non-positive, it's an invalid bet (leads to ruin)
        # Return a large value (infinity) to have the optimizer avoid this region.
        if g_A <= 1e-9 or g_B <= 1e-9 or g_C <= 1e-9:
            return np.inf

        # Calculate the expected log growth
        growth = p[0] * np.log(g_A) + p[1] * np.log(g_B) + p[2] * np.log(g_C)
        
        # Return the negative because we are using a minimizer
        return -growth

    # --- Step 1: Calculate the optimal fractions (f*) and growth rate (W*) ---
    
    # Constraints: f_i >= 0 (bounds) and sum(f_i) <= 1
    bounds = ((0, 1), (0, 1), (0, 1))
    constraints = ({'type': 'ineq', 'fun': lambda f: 1 - np.sum(f)})
    
    # Initial guess for the optimizer
    f_initial_guess = np.array([0.1, 0.1, 0.1])
    
    # Run the optimizer to find optimal fractions f* using TRUE probabilities
    opt_result_true = minimize(
        fun=neg_log_growth,
        x0=f_initial_guess,
        args=(p_true, b),
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        tol=1e-12
    )

    # The optimal growth rate W* is the negative of the minimized function value
    W_star = -opt_result_true.fun

    # --- Step 2: Calculate the fractions (f') based on mistaken beliefs ---
    
    # Run the optimizer to find fractions f' using MISTAKEN probabilities
    opt_result_believed = minimize(
        fun=neg_log_growth,
        x0=f_initial_guess,
        args=(p_believed, b),
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        tol=1e-12
    )
    
    # These are the fractions the person will actually bet
    f_actual_bets = opt_result_believed.x

    # --- Step 3: Calculate the actual growth rate (W) ---
    
    # The actual growth is calculated using the mistaken bets (f_actual_bets)
    # but with the outcomes governed by the true probabilities (p_true)
    W_actual = -neg_log_growth(f_actual_bets, p_true, b)

    # --- Step 4: Output the results ---
    
    difference = W_star - W_actual
    
    print("This equation will be solved: W* - W")
    print(f"Optimal growth rate (W*): {W_star:.6f}")
    print(f"Actual growth rate (W): {W_actual:.6f}")
    print(f"The difference (W* - W) is: {W_star:.6f} - {W_actual:.6f} = {difference:.6f}")
    
    return difference

# Execute the function and capture the final answer for submission format
final_answer = solve_race_growth()
print(f"<<<{final_answer:.6f}>>>")