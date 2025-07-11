import math

def calculate_payoffs():
    """
    Solves the game theory problem to find the value of 1-p.
    """
    N = 8
    n_players = 3 * N
    n_opponents = n_players - 1

    # Step 1: Define payoff for a player using Discrete Strategy (S_D)
    # against n_opponents also playing Discrete.
    # This is pi(S_D | K=23).
    pi_D_vs_23D = (1/3) * (1 - (1 - 1/N)**n_players)
    
    # Step 2: Define payoff for a player using Continuous Strategy S_m
    # against n_opponents playing Discrete.
    # The optimal deviation is found to be m=3. This is pi(S_3 | K=23).
    m = 3
    pi_C3_vs_23D = sum(
        ((-1)**(k - 1)) * math.comb(m, k) * (1 - k/N)**n_opponents
        for k in range(1, m + 1)
    )
    
    # Calculate Delta(23), the payoff difference in a mostly discrete world.
    # The premise is that pi_C3 > pi_D, so this is negative.
    delta_23 = pi_D_vs_23D - pi_C3_vs_23D

    # Step 3: Calculate payoffs against a mixed population of
    # 22 S_D players and 1 S_3 player (K=22, J=1).
    K = 22
    J = 1

    # Payoff for S_D against 22 S_D players. S_3 players are irrelevant
    # as their bid (1/3) is too low.
    pi_D_vs_22D_1C = (N / (K + 1)) * (1 - (1 - 1/N)**(K + 1))

    # Payoff for S_3 against 22 S_D and 1 S_3 opponent.
    # This requires an inclusion-exclusion calculation for winning any of the 3 races.
    
    # Prob of winning any single race (e.g., race 1)
    # This happens if no D-player chose it, and we win the tie-break against the C-player.
    prob_no_D_on_race1 = (1 - 1/N)**K
    # The C-opponent chooses 3 races. Prob they choose race 1 is 3/8.
    # Prob of winning the tie-break vs the single C-opponent for a single race:
    prob_win_tiebreak_C_1race = (1 - m/N) * 1 + (m/N) * (1/2) # (5/8)*1 + (3/8)*0.5 = 13/16
    P_W1 = prob_no_D_on_race1 * prob_win_tiebreak_C_1race
    
    # Prob of winning two specific races (e.g., race 1 & 2)
    prob_no_D_on_race12 = (1 - 2/N)**K
    # Prob of winning tie-breaks for 2 races against the C-opponent.
    # The C-opponent chooses 3 out of 8 races (56 ways).
    # Overlap with our 2 races can be 0, 1, or 2.
    # k=0 (C chooses 3 from other 6; 20 ways): tie-win prob=1*1=1
    # k=1 (C chooses 1 of our 2, 2 of other 6; 2*15=30 ways): tie-win prob=1*0.5=0.5
    # k=2 (C chooses 2 of our 2, 1 of other 6; 1*6=6 ways): tie-win prob=0.5*0.5=0.25
    # Since we assume simultaneous resolution, we need to win both.
    prob_C_no_overlap_2 = math.comb(N - 2, m) / math.comb(N, m) # 20/56
    prob_C_1_overlap_2 = (math.comb(2, 1) * math.comb(N-2, m-1)) / math.comb(N,m) # 30/56
    prob_C_2_overlap_2 = (math.comb(2, 2) * math.comb(N-2, m-2)) / math.comb(N,m) # 6/56
    prob_win_tiebreak_C_2races = (prob_C_no_overlap_2 * 1 
                                 + prob_C_1_overlap_2 * 0.5 
                                 + prob_C_2_overlap_2 * 0.25) # simplified from E[1/(k1+1) * 1/(k2+1)] logic
    P_W12 = prob_no_D_on_race12 * prob_win_tiebreak_C_2races

    # Prob of winning three specific races (race 1, 2, 3)
    prob_no_D_on_race123 = (1 - 3/N)**K
    # Prob of winning tie-breaks for 3 races against the C-opponent.
    # C-opponent chooses 3 races. They can overlap with ours by 0,1,2,3.
    prob_C_0_overlap_3 = math.comb(N-3,3)/math.comb(N,3) # 10/56
    prob_C_1_overlap_3 = (math.comb(3,1)*math.comb(N-3,2))/math.comb(N,3) # 30/56
    prob_C_2_overlap_3 = (math.comb(3,2)*math.comb(N-3,1))/math.comb(N,3) # 15/56
    prob_C_3_overlap_3 = (math.comb(3,3)*math.comb(N-3,0))/math.comb(N,3) # 1/56
    prob_win_tiebreak_C_3races = (prob_C_0_overlap_3 * 1 
                                 + prob_C_1_overlap_3 * 0.5 
                                 + prob_C_2_overlap_3 * 0.25
                                 + prob_C_3_overlap_3 * 0.125)
    P_W123 = prob_no_D_on_race123 * prob_win_tiebreak_C_3races
    
    pi_C3_vs_22D_1C = 3 * P_W1 - 3 * P_W12 + P_W123

    delta_22 = pi_D_vs_22D_1C - pi_C3_vs_22D_1C

    # Step 4: Calculate 1-p using the approximation
    one_minus_p = -delta_23 / (n_opponents * delta_22)
    
    # Step 5: Compute the final answer
    final_answer = math.floor(10000 * one_minus_p)
    
    print(f"For N = {N}, the optimal deviation strategy is m = {m}.")
    print(f"Payoff for Discrete player (vs all-Discrete): {pi_D_vs_23D:.8f}")
    print(f"Payoff for Continuous player (vs all-Discrete): {pi_C3_vs_23D:.8f}")
    print(f"The payoff difference, Delta(23), is {delta_23:.8f}.")
    print("-" * 20)
    print(f"Payoff for D (vs 22 D, 1 C): {pi_D_vs_22D_1C:.8f}")
    print(f"Payoff for C (vs 22 D, 1 C): {pi_C3_vs_22D_1C:.8f}")
    print(f"The payoff difference, Delta(22), is {delta_22:.8f}.")
    print("-" * 20)
    print(f"The calculated value for (1-p) is approximately {one_minus_p:.8f}.")
    print(f"The final calculation is floor(10000 * (1-p)) = floor({10000 * one_minus_p:.4f}).")
    print(f"The result is {final_answer}.")
    
    # According to the problem statement, we should output the final equation with all numbers.
    print("\nFinal Equation:")
    print(f"p is the solution to the equilibrium condition. We approximate 1-p.")
    print(f"1-p ≈ - (Δ(23) / (23 * Δ(22)))")
    print(f"1-p ≈ - (({pi_D_vs_23D:.6f} - {pi_C3_vs_23D:.6f}) / (23 * ({pi_D_vs_22D_1C:.6f} - {pi_C3_vs_22D_1C:.6f})))")
    print(f"1-p ≈ - ({delta_23:.6f} / (23 * {delta_22:.6f}))")
    print(f"1-p ≈ {one_minus_p:.6f}")
    print(f"floor(10000 * {one_minus_p:.6f}) = {final_answer}")


calculate_payoffs()