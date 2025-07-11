import collections

def solve_for_p():
    """
    This function calculates the value of p based on the provided economic model and data.
    """
    # Step 1: Define the given probability distribution and deduce the states of the world.
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # From y = s * e, we can find the values of s.
    # For e=22, y=132 -> s = 132 / 22 = 6
    # For e=22, y=44  -> s = 44 / 22 = 2
    s_H = 6
    s_L = 2
    print(f"Step 1: From the data, we deduce the states of the world: s_H = {s_H}, s_L = {s_L}")

    # Step 2: Map the (e, y) distribution to an (e, s) distribution.
    # We associate the higher effort e=22 with a 'high' signal and e=10 with a 'low' signal.
    P_e22_sH = prob_dist[(22, 132)]
    P_e22_sL = prob_dist[(22, 44)]
    P_e10_sH = prob_dist[(10, 60)]
    P_e10_sL = prob_dist[(10, 20)]

    # Step 3: Calculate conditional expectations of s given the effort choice (signal).
    # First, find the marginal probability of each effort level.
    P_e22 = P_e22_sH + P_e22_sL
    P_e10 = P_e10_sH + P_e10_sL

    # Then, find the conditional probability of the state given the effort.
    P_sH_given_e22 = P_e22_sH / P_e22
    P_sL_given_e22 = P_e22_sL / P_e22
    P_sH_given_e10 = P_e10_sH / P_e10
    P_sL_given_e10 = P_e10_sL / P_e10

    # Now, compute the conditional expectations.
    E_s_given_e22 = s_H * P_sH_given_e22 + s_L * P_sL_given_e22
    E_s_given_e10 = s_H * P_sH_given_e10 + s_L * P_sL_given_e10
    print(f"\nStep 2: Calculate the employee's conditional expectations.")
    print(f"E[s | signal for e=22] = {s_H}*{P_sH_given_e22:.3f} + {s_L}*{P_sL_given_e22:.3f} = {E_s_given_e22}")
    print(f"E[s | signal for e=10] = {s_H}*{P_sH_given_e10:.3f} + {s_L}*{P_sL_given_e10:.3f} = {E_s_given_e10}")

    # Step 4: Determine β using the employee's optimal effort rule: e = β * E[s|signal].
    # We can use either effort level. Let's use e=22.
    e_high = 22
    beta = e_high / E_s_given_e22
    print(f"\nStep 3: Determine the contract parameter β.")
    print(f"From the employee's effort rule, we have the equation: {e_high} = β * {E_s_given_e22}")
    print(f"Solving for β gives: β = {e_high} / {E_s_given_e22} = {beta}")

    # Step 5 & 6: The firm's profit maximization implies p = β.
    # The first-order condition for the firm's choice of β simplifies to p = β.
    p = beta
    print("\nStep 4: Determine the price p.")
    print("The firm chooses β to maximize profit. The first-order condition of this problem simplifies to p = β.")
    print(f"\nFinal Equation: p = β")
    print(f"The value of p is {p}.")

solve_for_p()