import collections

def solve_price():
    """
    Calculates the price 'p' based on the provided probability distribution of effort and output.
    """
    # Step 1: Define the given probability distribution P(e, y)
    # The data is given as P(e, y) pairs.
    # We can map this to P(e, s) using s = y/e.
    # e_H = 22, e_L = 10
    # s_H = 6, s_L = 2
    prob_e_s = {
        (22, 6): 0.4375,  # P(e=22, y=132)
        (22, 2): 0.0625,  # P(e=22, y=44)
        (10, 6): 0.0625,  # P(e=10, y=60)
        (10, 2): 0.4375,  # P(e=10, y=20)
    }

    e_H = 22
    e_L = 10
    s_H = 6
    s_L = 2

    # Step 2: Calculate marginal probabilities for effort levels
    prob_e_H = prob_e_s[(e_H, s_H)] + prob_e_s[(e_H, s_L)]
    prob_e_L = prob_e_s[(e_L, s_H)] + prob_e_s[(e_L, s_L)]

    # Step 3: Calculate the conditional expectation of the state 's' given high effort e_H
    # E[s | e=e_H] = sum(s * P(s | e=e_H)) = sum(s * P(s, e=e_H) / P(e=e_H))
    expected_s_given_e_H_num = s_H * prob_e_s[(e_H, s_H)] + s_L * prob_e_s[(e_H, s_L)]
    expected_s_given_e_H = expected_s_given_e_H_num / prob_e_H
    
    # Step 4: Calculate beta using the employee's first-order condition: e_H = beta * E[s|e=e_H]
    beta = e_H / expected_s_given_e_H

    # Step 5: Based on principal-agent theory, for a risk-neutral agent, the optimal contract
    # sets the incentive `beta` equal to the price `p` to maximize the firm's profit.
    p = beta

    # Print the explanation and the final equation with numbers
    print("Step 1: Calculate the employee's expected state of the world given high effort (e=22).")
    print(f"E[s|e={e_H}] = (s_H * P(s={s_H},e={e_H}) + s_L * P(s={s_L},e={e_H})) / P(e={e_H})")
    print(f"E[s|e={e_H}] = ({s_H} * {prob_e_s[(e_H, s_H)]} + {s_L} * {prob_e_s[(e_H, s_L)]}) / {prob_e_H}")
    print(f"E[s|e={e_H}] = {expected_s_given_e_H_num} / {prob_e_H} = {expected_s_given_e_H}")
    print("\nStep 2: Calculate the contract parameter 'beta' from the employee's effort choice.")
    print("The employee's optimal effort rule is e = beta * E[s|signal].")
    print(f"beta = e_H / E[s|e={e_H}]")
    print(f"beta = {e_H} / {expected_s_given_e_H} = {beta}")
    print("\nStep 3: Determine the price 'p'.")
    print("In this principal-agent model, the firm's profit-maximizing strategy is to set p = beta.")
    print("\nThe final equation is p = e_H / E[s|e=e_H]")
    print(f"{p} = {e_H} / {expected_s_given_e_H}")
    print(f"\nThe value of p is {p}.")

solve_price()