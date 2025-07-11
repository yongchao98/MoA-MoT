import sys

def solve():
    """
    Solves for the price p based on the given economic model and data.
    """
    # The given joint probability distribution of (effort, output) pairs
    P_e_y = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375,
    }

    # High effort e_H corresponds to the employee observing a high signal (θ=s_H)
    e_H = 22
    # Low effort e_L corresponds to the employee observing a low signal (θ=s_L)
    e_L = 10
    
    # Step 1: Determine the states of the world, s_H and s_L
    # From y = s * e, we can find the state s = y / e
    # When e=22, a high output (132) implies a high state (s_H)
    # and a low output (44) implies a low state (s_L).
    s_H = int(132 / e_H)
    s_L = int(44 / e_H)

    # We can check for consistency using the other effort level
    # s_H check: 60 / e_L = 6
    # s_L check: 20 / e_L = 2
    
    print("Step 1: Determine the states of the world (s_H, s_L).")
    print(f"The high state is s_H = 132 / {e_H} = {s_H}")
    print(f"The low state is s_L = 44 / {e_H} = {s_L}")
    print("-" * 30)

    # Step 2: Calculate the probability of observing each signal
    # P(θ=s_H) is the marginal probability of the event e=22
    P_theta_H = P_e_y[(22, 132)] + P_e_y[(22, 44)]
    # P(θ=s_L) is the marginal probability of the event e=10
    P_theta_L = P_e_y[(10, 60)] + P_e_y[(10, 20)]

    print("Step 2: Calculate signal probabilities P(θ).")
    print(f"P(θ=s_H) = P(e={e_H}) = {P_e_y[(22, 132)]} + {P_e_y[(22, 44)]} = {P_theta_H}")
    print(f"P(θ=s_L) = P(e={e_L}) = {P_e_y[(10, 60)]} + {P_e_y[(10, 20)]} = {P_theta_L}")
    print("-" * 30)

    # Step 3: Calculate conditional probabilities P(s|θ) using Bayes' rule
    # P(s=s_H | θ=s_H) = P(s=s_H, θ=s_H) / P(θ=s_H)
    # The event (s=s_H, θ=s_H) corresponds to (e=22, y=132)
    P_sH_given_thetaH = P_e_y[(22, 132)] / P_theta_H
    P_sL_given_thetaH = P_e_y[(22, 44)] / P_theta_H
    P_sH_given_thetaL = P_e_y[(10, 60)] / P_theta_L
    P_sL_given_thetaL = P_e_y[(10, 20)] / P_theta_L
    
    print("Step 3: Calculate conditional probabilities P(s|θ).")
    print(f"P(s=s_H | θ=s_H) = {P_e_y[(22, 132)]} / {P_theta_H} = {P_sH_given_thetaH}")
    print(f"P(s=s_L | θ=s_H) = {P_e_y[(22, 44)]} / {P_theta_H} = {P_sL_given_thetaH}")
    print(f"P(s=s_H | θ=s_L) = {P_e_y[(10, 60)]} / {P_theta_L} = {P_sH_given_thetaL}")
    print(f"P(s=s_L | θ=s_L) = {P_e_y[(10, 20)]} / {P_theta_L} = {P_sL_given_thetaL}")
    print("-" * 30)

    # Step 4: Calculate conditional expectations E[s|θ]
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

    print("Step 4: Calculate conditional expected states E[s|θ].")
    print(f"E[s | θ=s_H] = {s_H}*{P_sH_given_thetaH} + {s_L}*{P_sL_given_thetaH} = {E_s_given_thetaH}")
    print(f"E[s | θ=s_L] = {s_H}*{P_sH_given_thetaL} + {s_L}*{P_sL_given_thetaL} = {E_s_given_thetaL}")
    print("-" * 30)

    # Step 5: Solve for beta using the employee's effort rule: e = β * E[s|θ]
    beta = int(e_H / E_s_given_thetaH)
    
    print("Step 5: Determine the contract parameter β.")
    print(f"From e_H = β * E[s|θ=s_H] => β = {e_H} / {E_s_given_thetaH} = {beta}")
    # Consistency Check
    print(f"Check with e_L = β * E[s|θ=s_L] => β = {e_L} / {E_s_given_thetaL} = {int(e_L / E_s_given_thetaL)}")
    print("-" * 30)

    # Step 6: Use the firm's profit maximization condition to find p.
    # The firm chooses β to maximize E[Profit] = p*E[y] - E[e^2/2].
    # The first-order condition for this maximization simplifies to p = β.
    p = beta
    
    print("Step 6: Determine the price p.")
    print("The firm's optimal choice of β satisfies the condition p = β.")
    print(f"Since β = {beta}, p must be equal to {p}.")
    print("\nFinal equation:")
    print(f"{p} = {beta}")

solve()