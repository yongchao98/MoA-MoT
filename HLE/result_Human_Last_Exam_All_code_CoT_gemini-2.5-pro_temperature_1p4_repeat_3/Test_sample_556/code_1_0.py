import sys

# Suppress writing bytecode files
sys.dont_write_bytecode = True

def solve_price():
    """
    This function solves for the price p based on the provided economic model and data.
    The steps are:
    1. Determine the states of the world (s_H, s_L).
    2. Calculate the conditional expectations of the state, E[s | θ].
    3. Solve for the incentive parameter β using the employee's optimal effort rule.
    4. Equate p to β based on the firm's profit maximization problem.
    """

    # Given probabilities from the problem description
    P_e22_y132 = 0.4375  # This is P(effort=high, state=high)
    P_e22_y44 = 0.0625   # This is P(effort=high, state=low)
    P_e10_y60 = 0.0625   # This is P(effort=low, state=high)
    P_e10_y20 = 0.4375   # This is P(effort=low, state=low)

    # Observed effort levels and corresponding outputs
    e_H = 22.0  # Effort when signal is high
    e_L = 10.0  # Effort when signal is low
    y_H_H = 132.0 # Output when effort is high and state is high
    y_H_L = 44.0  # Output when effort is high and state is low
    
    print("Step 1: Determine the states of the world (s_H and s_L).")
    print("The output y is defined as y = s * e, where s is the state and e is the effort.")
    
    # Calculate s_H and s_L from the data where effort is high (e_H)
    s_H = y_H_H / e_H
    s_L = y_H_L / e_H
    
    print(f"Using the data for high effort (e={e_H}):")
    print(f"The high state s_H = output / effort = {y_H_H} / {e_H} = {s_H}")
    print(f"The low state s_L = output / effort = {y_H_L} / {e_H} = {s_L}")
    print("-" * 40)

    print("Step 2: Calculate the employee's conditional expectation of the state E[s|θ].")
    print("The employee's effort e depends on the signal θ. We can proxy the signal with the observed effort level.")
    
    # Calculate the marginal probability of observing a high or low effort
    # P(e=high) is equivalent to P(θ=s_H)
    P_e_high = P_e22_y132 + P_e22_y44
    # P(e=low) is equivalent to P(θ=s_L)
    P_e_low = P_e10_y60 + P_e10_y20
    
    print(f"The probability of observing high effort is P(e={e_H}) = {P_e22_y132} + {P_e22_y44} = {P_e_high}")
    print(f"The probability of observing low effort is P(e={e_L}) = {P_e10_y60} + {P_e10_y20} = {P_e_low}")

    # Calculate conditional probabilities P(s|e) using P(s,e) / P(e)
    # When effort is high
    P_sH_given_eH = P_e22_y132 / P_e_high
    P_sL_given_eH = P_e22_y44 / P_e_high
    
    # When effort is low
    P_sH_given_eL = P_e10_y60 / P_e_low
    P_sL_given_eL = P_e10_y20 / P_e_low
    
    # Calculate conditional expectations E[s|e]
    E_s_given_eH = s_H * P_sH_given_eH + s_L * P_sL_given_eH
    E_s_given_eL = s_H * P_sH_given_eL + s_L * P_sL_given_eL

    print(f"The conditional expectation of s given high effort is E[s|e={e_H}] = {s_H}*{P_sH_given_eH:.3f} + {s_L}*{P_sL_given_eH:.3f} = {E_s_given_eH}")
    print(f"The conditional expectation of s given low effort is E[s|e={e_L}] = {s_H}*{P_sH_given_eL:.3f} + {s_L}*{P_sL_given_eL:.3f} = {E_s_given_eL}")
    print("-" * 40)

    print("Step 3: Solve for the incentive parameter β.")
    print("The employee chooses effort e to maximize E[w - e^2/2 | θ], which gives the optimal effort rule: e = β * E[s|θ].")
    
    # Calculate β using the high effort case
    beta = e_H / E_s_given_eH
    # Verify with the low effort case
    beta_check = e_L / E_s_given_eL

    print(f"Using the high effort case: {e_H} = β * {E_s_given_eH}")
    print(f"This gives β = {e_H} / {E_s_given_eH} = {beta}")
    print(f"Checking with the low effort case: {e_L} = β * {E_s_given_eL} gives β = {beta_check}, confirming the value.")
    print("-" * 40)

    print("Step 4: Solve for the price p.")
    print("The firm chooses β to maximize its expected profit. In this principal-agent model, the profit-maximizing choice is β = p.")
    
    # The price p must be equal to the optimal beta
    p = beta
    
    print(f"Since the observed data results from this optimization, the price p must be equal to β.")
    print(f"Therefore, p = {p}")
    print("-" * 40)

    print("The final equation can be written in terms of the optimal effort rule, substituting p for β:")
    print(f"e_H = p * E[s|e_H]")
    print(f"{int(e_H)} = {int(p)} * {E_s_given_eH}")

    return p

if __name__ == '__main__':
    price = solve_price()
    print(f"<<<{price}>>>")