import collections

def solve_for_p():
    """
    Calculates the price p based on the given probability distribution of effort and output.
    """
    # Step 1: Define the given data and identify model parameters.
    # P(e, y) is the joint probability distribution.
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # The effort levels are 22 and 10. Let's assume e_H > e_L.
    e_H = 22
    e_L = 10

    # From y = s * e, we can find the states of the world, s.
    # For e = 22, y = 132 -> s = 132 / 22 = 6
    # For e = 22, y = 44  -> s = 44 / 22 = 2
    # So, s_H = 6 and s_L = 2.

    s_H = 6
    s_L = 2

    # Step 2: Calculate the conditional expectation E[s | e].
    # We will use the high effort case, e_H = 22. The signal associated with this
    # effort is θ_H. We need to calculate E[s | θ_H], which is equivalent to E[s | e=22].
    
    # Probability of events given e=22
    prob_eH_sH = prob_dist[(e_H, e_H * s_H)] # P(e=22, s=6)
    prob_eH_sL = prob_dist[(e_H, e_H * s_L)] # P(e=22, s=2)

    # Marginal probability of choosing high effort e_H
    prob_eH = prob_eH_sH + prob_eH_sL

    # Conditional probabilities P(s | e=22)
    prob_sH_given_eH = prob_eH_sH / prob_eH
    prob_sL_given_eH = prob_eH_sL / prob_eH

    # Conditional expectation E[s | e=22]
    E_s_given_eH = s_H * prob_sH_given_eH + s_L * prob_sL_given_eH
    
    # Step 3: Solve for p using the relationship e = p * E[s | e].
    # For the high effort case: e_H = p * E[s | e=e_H]
    p = e_H / E_s_given_eH

    # Print the calculation steps and the final answer.
    print(f"The high effort level is e_H = {e_H}.")
    print(f"The conditional expectation of the state given high effort is E[s|e_H] = {E_s_given_eH}.")
    print(f"The price p is determined by the equation: e_H = p * E[s|e_H].")
    print(f"Therefore, p = {e_H} / {E_s_given_eH} = {p}")

solve_for_p()