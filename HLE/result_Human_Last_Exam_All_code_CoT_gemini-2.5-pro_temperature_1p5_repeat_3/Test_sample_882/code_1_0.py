import math

def solve_game_theory_problem():
    """
    This function solves the game theory problem to find the mixed strategy
    equilibrium probability 'p' and then computes the final requested value.
    """

    # 1. Define constants for N=8
    N = 8
    M = 3 * N - 1  # Number of other players
    X = 1.0 - 1.0 / N # Probability of a discrete player NOT choosing a specific race

    # 2. Define the function f(p) = pi_D(p) - pi_S(p) to find its root.
    # pi_D(p) is the expected payoff (wins) for the Discrete strategy.
    # pi_S(p) is the expected payoff (wins) for the Spread strategy.
    # The function calculates the expectation over the binomial distribution of
    # d (number of discrete opponents).
    def f(p):
        if p < 0 or p > 1:
            return float('nan')
        # Handle edge cases for numerical stability and correctness
        if p == 0:
            # Only d=0 term is non-zero
            payoff_d_discrete = (N / 1.0) * (1.0 - X**1)
            payoff_s_spread = N * (X**0) / (M + 1.0)
            return payoff_d_discrete - payoff_s_spread
        if p == 1:
            # Only d=M term is non-zero
            payoff_d_discrete = (N / (M + 1.0)) * (1.0 - X**(M + 1))
            payoff_s_spread = N * (X**M) / 1.0
            return payoff_d_discrete - payoff_s_spread

        total_diff = 0.0
        log_p = math.log(p)
        log_1_minus_p = math.log(1.0 - p)

        for d in range(M + 1):
            # Binomial probability P(d) using log-gamma for precision
            log_binom_coeff = math.lgamma(M + 1) - math.lgamma(d + 1) - math.lgamma(M - d + 1)
            log_prob_d = log_binom_coeff + d * log_p + (M - d) * log_1_minus_p
            prob_d = math.exp(log_prob_d)

            # Payoffs given d other players chose discrete
            # Payoff for player choosing Discrete strategy
            payoff_d_discrete = (N / (d + 1.0)) * (1.0 - X**(d + 1))
            
            # Payoff for player choosing Spread strategy
            payoff_s_spread = N * (X**d) / (M - d + 1.0)
            
            diff = payoff_d_discrete - payoff_s_spread
            total_diff += prob_d * diff
            
        return total_diff

    # 3. Use Bisection method to find the root p
    low_p = 0.0
    high_p = 1.0
    
    # We need high precision for p to ensure the final integer is correct.
    # 100 iterations are more than sufficient.
    for _ in range(100):
        mid_p = (low_p + high_p) / 2.0
        f_mid = f(mid_p)
        
        # f(p) is a decreasing function, f(0) > 0, f(1) < 0
        if f_mid > 0:
            low_p = mid_p
        else:
            high_p = mid_p

    p_equilibrium = (low_p + high_p) / 2.0

    # 4. Calculate and print the final result as requested
    one_minus_p = 1.0 - p_equilibrium
    value_to_floor = 10000 * one_minus_p
    final_answer = math.floor(value_to_floor)
    
    # Format for printing, using 6 significant digits for the value 1-p
    p_val_str = f"{p_equilibrium:.6f}"
    one_minus_p_str = f"{one_minus_p:.6f}" 
    value_to_floor_str = f"{value_to_floor:.6f}"

    print(f"The equilibrium probability 'p' that a player devotes all fuel to a single race is found by numerically solving the indifference equation.\n")
    print(f"The solution for N=8 is p = {p_equilibrium:.8f}\n")
    print(f"The problem asks for the value of floor(10000 * (1 - p)).")
    print(f"We will now compute the value of each term in this equation:\n")
    print(f"10000")
    print(f"1")
    print(f"p = {p_val_str}\n")
    print(f"The step-by-step calculation is as follows:")
    print(f"floor(10000 * (1 - {p_val_str}))")
    print(f"= floor(10000 * {one_minus_p_str})")
    print(f"= floor({value_to_floor_str})")
    print(f"= {final_answer}")
    
    return final_answer

# Execute the solver and print the final answer in the required format
result = solve_game_theory_problem()
print(f"<<<{result}>>>")
