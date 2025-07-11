import sys

def solve():
    """
    This script calculates the price 'p' based on the provided probability distribution
    of an employee's effort and output in a principal-agent model.
    """
    
    # 1. Define known values from the problem statement
    # P(e, y) pairs are interpreted as P(signal, state)
    P_thetaH_sH = 0.4375  # P(e=22, y=132) -> P(θ=s_H, s=s_H)
    P_thetaH_sL = 0.0625  # P(e=22, y=44)  -> P(θ=s_H, s=s_L)
    P_thetaL_sH = 0.0625  # P(e=10, y=60)  -> P(θ=s_L, s=s_H)
    P_thetaL_sL = 0.4375  # P(e=10, y=20)  -> P(θ=s_L, s=s_L)

    e_H = 22  # High effort, corresponding to signal θ=s_H
    e_L = 10  # Low effort, corresponding to signal θ=s_L
    
    y_H_H = 132 # Output when e=e_H and s=s_H

    # 2. Determine the states of the world, s_H and s_L
    # Since y = s * e, we can find s. We assume s_H > s_L.
    # When effort is high (e_H=22), output can be 132 or 44.
    s_H = 132 / 22
    s_L = 44 / 22

    print("Step 1: Determine the states of the world (s_H, s_L).")
    print(f"From the data where e = 22, the possible states are s = 132/22 = {s_H} and s = 44/22 = {s_L}.")
    print(f"So, s_H = {s_H} and s_L = {s_L}.")
    print("-" * 30)

    # 3. Calculate conditional expectations E[s|θ]
    # First, calculate marginal probabilities of the signals P(θ)
    P_theta_H = P_thetaH_sH + P_thetaH_sL
    P_theta_L = P_thetaL_sH + P_thetaL_sL

    # Next, calculate conditional probabilities of the state P(s|θ) using Bayes' rule.
    # For θ = s_H:
    P_sH_given_thetaH = P_thetaH_sH / P_theta_H
    P_sL_given_thetaH = P_thetaH_sL / P_theta_H
    E_s_given_theta_H = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH

    # For θ = s_L:
    P_sH_given_thetaL = P_thetaL_sH / P_theta_L
    P_sL_given_thetaL = P_thetaL_sL / P_theta_L
    E_s_given_theta_L = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL
    
    print("Step 2: Calculate the worker's conditional expectations of the state.")
    print(f"The worker's expected state given a high signal is E[s|θ=s_H] = {s_H}*({P_thetaH_sH}/{P_theta_H}) + {s_L}*({P_thetaH_sL}/{P_theta_H}) = {E_s_given_theta_H}")
    print(f"The worker's expected state given a low signal is E[s|θ=s_L] = {s_H}*({P_thetaL_sH}/{P_theta_L}) + {s_L}*({P_thetaL_sL}/{P_theta_L}) = {E_s_given_theta_L}")
    print("-" * 30)

    # 4. Determine the contract parameter β
    # The employee's optimal effort is e = β * E[s|θ].
    beta_from_H = e_H / E_s_given_theta_H
    beta_from_L = e_L / E_s_given_theta_L
    
    # The value of beta should be consistent
    if not (0.999 * beta_from_L < beta_from_H < 1.001 * beta_from_L):
        print("Error: Inconsistent values of beta derived. Check logic.", file=sys.stderr)
        return

    beta = beta_from_H

    print("Step 3: Determine the contract parameter β.")
    print(f"From the incentive compatibility condition e = β * E[s|θ]:")
    print(f"Using the high effort case: {e_H} = β * {E_s_given_theta_H:.1f}, which gives β = {beta_from_H:.1f}")
    print(f"Using the low effort case: {e_L} = β * {E_s_given_theta_L:.1f}, which gives β = {beta_from_L:.1f}")
    print(f"The value is consistent. β = {beta}.")
    print("-" * 30)

    # 5. Determine the price p
    # The firm chooses β to maximize profit. The first-order condition for this model simplifies to p = β.
    # Therefore, the market price p must equal the optimal β chosen by the firm.
    p = beta

    print("Step 4: Determine the price p based on the firm's profit maximization.")
    print("The firm's optimal choice of the contract parameter β must satisfy the condition that maximizes profit.")
    print("For a linear contract and quadratic effort cost, the maximization of total surplus leads to the condition: p = β.")
    print("The final equation is p = β.")
    print(f"So, p = {p:.0f}")

    # Final Answer
    # The following line is for the final answer format.
    print(f"\n<<<{p:.0f}>>>")

solve()