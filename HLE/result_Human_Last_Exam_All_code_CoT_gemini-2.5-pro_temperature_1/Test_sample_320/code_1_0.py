from fractions import Fraction

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed of a 3-particle system.
    """
    # Step 1: Set up the system of equations for q1 = P(y1>1) and q2 = P(y2>1)
    # from the condition that all particle velocities must be equal in steady state.
    # v1 = q1 - 1/3
    # v2 = q2 - q1
    # v3 = 1 - q2
    #
    # Setting v1 = v3 gives: q1 - 1/3 = 1 - q2  => q1 + q2 = 4/3  (Eq. 1)
    # Setting v1 = v2 gives: q1 - 1/3 = q2 - q1  => 2*q1 - q2 = 1/3 (Eq. 2)
    
    # We solve this 2x2 linear system.
    # Adding (Eq. 1) and (Eq. 2):
    # (q1 + q2) + (2*q1 - q2) = 4/3 + 1/3
    # 3*q1 = 5/3  => q1 = 5/9
    
    q1_val = Fraction(5, 9)
    
    # Substitute q1 into (Eq. 1):
    # 5/9 + q2 = 4/3 => q2 = 4/3 - 5/9 = 12/9 - 5/9 = 7/9
    q2_val = Fraction(7, 9)

    print("--- Calculating Probabilities and Speed ---")
    print(f"The probability of the first gap being greater than 1 is q1 = P(y1 > 1) = {q1_val}")
    print(f"The probability of the second gap being greater than 1 is q2 = P(y2 > 1) = {q2_val}")
    
    # Step 2: Calculate the asymptotic speed v using any of the velocity equations.
    # v = q1 - 1/3
    rate_leftmost_L = Fraction(1, 3)
    speed = q1_val - rate_leftmost_L
    
    print(f"\nThe asymptotic speed of the leftmost particle is v = q1 - 1/3 = {q1_val} - {rate_leftmost_L} = {speed}")

    # Step 3: Calculate average distances.
    # The stationary distribution for the gaps (y1, y2) is a product of two
    # geometric distributions with parameters alpha = q1 and beta = q2.
    # P(y_i=k) = (1-p) * p^(k-1) for k>=1.
    # The mean of such a distribution is 1 / (1-p).
    
    alpha = q1_val
    beta = q2_val
    
    # E[y1] = 1 / (1 - alpha)
    E_y1 = 1 / (1 - alpha)
    
    # E[y2] = 1 / (1 - beta)
    E_y2 = 1 / (1 - beta)
    
    print("\n--- Calculating Average Distances ---")
    print(f"The distribution for the first gap y1 is geometric with parameter alpha = {alpha}.")
    print(f"The average first gap is E[y1] = 1 / (1 - alpha) = 1 / (1 - {alpha}) = {E_y1}")
    
    print(f"\nThe distribution for the second gap y2 is geometric with parameter beta = {beta}.")
    print(f"The average second gap is E[y2] = 1 / (1 - beta) = 1 / (1 - {beta}) = {E_y2}")

    # The total average distance is the sum of the average gaps.
    total_distance = E_y1 + E_y2
    
    print(f"\nThe average distance between the leftmost and rightmost particles is E[y1] + E[y2] = {E_y1} + {E_y2} = {total_distance}")
    
    # Step 4: Final Answer
    print("\n--- Final Answer ---")
    final_answer = (total_distance, speed)
    print(f"The result (distance, speed) is {final_answer}")
    
solve_particle_problem()