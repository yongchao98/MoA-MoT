import numpy as np
from scipy.optimize import minimize, LinearConstraint

def solve_three_horse_race():
    """
    Solves the Kelly Criterion problem for a three-horse race.

    The problem asks for W* - W, where:
    - W* is the optimal growth rate using true probabilities.
    - W is the actual growth rate from using a mistaken set of probabilities.
    """

    # --- Problem Setup ---
    # True probabilities for competitors (A, B, C)
    p_true = np.array([1/2, 1/4, 1/4])
    # Mistakenly believed probabilities
    p_mistaken = np.array([1/4, 1/2, 1/4])
    # Payout ratios (net odds b) for A, B, C
    b = np.array([4, 3, 3])

    # --- Growth Rate Calculation Function ---
    def calculate_growth_rate(fractions, probs, odds):
        """
        Calculates the expected log growth rate G = E[ln(Wealth_1 / Wealth_0)].
        
        Args:
            fractions (np.array): The fraction of capital bet on each competitor.
            probs (np.array): The probabilities of each competitor winning.
            odds (np.array): The net odds for each competitor.
            
        Returns:
            The growth rate, or a very small number if any wealth factor is non-positive.
        """
        f = np.array(fractions)
        F_total = np.sum(f)

        # Calculate the wealth factor R_i for each outcome i
        # R_i = 1 + f_i*b_i - sum_{j!=i} f_j
        # This can be rewritten as: R_i = 1 - F_total + f_i*(b_i + 1)
        wealth_factors = 1 - F_total + f * (odds + 1)
        
        # Log is undefined for non-positive wealth factors.
        if np.any(wealth_factors <= 0):
            return -np.inf # Return a very small number for invalid bets

        # Expected log growth rate
        g = np.sum(probs * np.log(wealth_factors))
        return g

    # --- Optimization Function ---
    def find_optimal_fractions(probs, odds):
        """
        Finds the optimal betting fractions f that maximize the growth rate.
        
        Args:
            probs (np.array): The probabilities used for optimization.
            odds (np.array): The net odds.
            
        Returns:
            The optimal fraction vector f.
        """
        # We want to maximize calculate_growth_rate, which is equivalent to
        # minimizing its negative.
        objective_func = lambda f: -calculate_growth_rate(f, probs, odds)
        
        # Number of competitors
        n_outcomes = len(probs)
        
        # Initial guess for the fractions
        f0 = np.ones(n_outcomes) / n_outcomes
        
        # Constraints:
        # 1. Each fraction must be non-negative: f_i >= 0
        #    This is handled by bounds.
        # 2. The sum of fractions cannot exceed 1: sum(f_i) <= 1
        bounds = [(0, 1)] * n_outcomes
        constraint_matrix = [[1] * n_outcomes]
        sum_constraint = LinearConstraint(constraint_matrix, lb=0, ub=1)

        result = minimize(
            objective_func,
            f0,
            method='SLSQP',
            bounds=bounds,
            constraints=[sum_constraint]
        )
        
        return result.x

    # --- Step 1: Calculate optimal strategy and growth rate (W*) ---
    f_optimal = find_optimal_fractions(p_true, b)
    W_star = calculate_growth_rate(f_optimal, p_true, b)
    
    # Calculate wealth factors for optimal strategy
    F_total_opt = np.sum(f_optimal)
    R_optimal = 1 - F_total_opt + f_optimal * (b + 1)
    
    # --- Step 2: Calculate sub-optimal strategy and actual growth rate (W) ---
    f_actual = find_optimal_fractions(p_mistaken, b)
    W_actual = calculate_growth_rate(f_actual, p_true, b)
    
    # Calculate wealth factors for actual (mistaken) strategy
    F_total_act = np.sum(f_actual)
    R_actual = 1 - F_total_act + f_actual * (b + 1)

    # --- Step 3: Display results ---
    print("This problem requires solving a constrained optimization problem to find the optimal betting fractions according to the Kelly criterion.")
    print("The goal is to maximize the expected logarithm of the wealth growth rate, G = p_A*ln(R_A) + p_B*ln(R_B) + p_C*ln(R_C), where R_i is the wealth factor if competitor i wins.\n")
    
    print("1. Optimal Strategy (based on true probabilities P=[0.500, 0.250, 0.250]):")
    print(f"   - Optimal betting fractions (f*): A={f_optimal[0]:.4f}, B={f_optimal[1]:.4f}, C={f_optimal[2]:.4f}")
    print(f"   - Optimal wealth factors (R*): R_A*={R_optimal[0]:.4f}, R_B*={R_optimal[1]:.4f}, R_C*={R_optimal[2]:.4f}")
    print( "   - Optimal growth rate equation:")
    print(f"     W* = {p_true[0]:.3f} * ln({R_optimal[0]:.4f}) + {p_true[1]:.3f} * ln({R_optimal[1]:.4f}) + {p_true[2]:.3f} * ln({R_optimal[2]:.4f})")
    print(f"     W* = {W_star:.6f}\n")

    print("2. Actual Strategy (based on mistaken probabilities P'=[0.250, 0.500, 0.250]):")
    print(f"   - Actual betting fractions (f): A={f_actual[0]:.4f}, B={f_actual[1]:.4f}, C={f_actual[2]:.4f}")
    print(f"   - Resulting wealth factors (R): R_A={R_actual[0]:.4f}, R_B={R_actual[1]:.4f}, R_C={R_actual[2]:.4f}")
    print( "   - Actual growth rate equation (using true probabilities with mistaken fractions):")
    print(f"     W = {p_true[0]:.3f} * ln({R_actual[0]:.4f}) + {p_true[1]:.3f} * ln({R_actual[1]:.4f}) + {p_true[2]:.3f} * ln({R_actual[2]:.4f})")
    print(f"     W = {W_actual:.6f}\n")
    
    difference = W_star - W_actual
    print("3. Difference (Optimal vs. Actual Growth Rate):")
    print(f"   W* - W = {W_star:.6f} - {W_actual:.6f} = {difference:.6f}")
    
    return difference

# Execute the solver and capture the final numerical answer
final_answer = solve_three_horse_race()