import collections

def solve_for_p():
    """
    This function calculates the price 'p' based on the provided probability distribution
    of effort and output in a principal-agent model.
    """
    # Step 1: Define the given probability distribution P(e, y)
    prob_dist = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    # Step 2: Determine the states of the world (s_H, s_L) and effort levels (e_H, e_L)
    # The production function is y = s * e. We can find 's' from the (e, y) pairs.
    # From (e=22, y=132), s = 132 / 22 = 6
    # From (e=22, y=44),  s = 44 / 22 = 2
    states = set()
    efforts = set()
    for (e, y), prob in prob_dist.items():
        states.add(y / e)
        efforts.add(e)

    s_H = max(states)
    s_L = min(states)
    e_H = max(efforts)
    e_L = min(efforts)

    print(f"From the data, we derive the model parameters:")
    print(f"High state of the world: s_H = {s_H}")
    print(f"Low state of the world: s_L = {s_L}")
    print(f"High effort level: e_H = {e_H}")
    print(f"Low effort level: e_L = {e_L}\n")

    # Step 3: Calculate conditional expectations E[s|e]
    # First, we need the joint probability P(e, s) and marginal probability P(e)
    # P(e, s) is derived from P(e, y) since y = s * e
    prob_e_s = {
        (e_H, s_H): prob_dist.get((e_H, e_H * s_H), 0),
        (e_H, s_L): prob_dist.get((e_H, e_H * s_L), 0),
        (e_L, s_H): prob_dist.get((e_L, e_L * s_H), 0),
        (e_L, s_L): prob_dist.get((e_L, e_L * s_L), 0)
    }

    # Marginal probability P(e)
    P_e_H = prob_e_s[(e_H, s_H)] + prob_e_s[(e_H, s_L)]
    P_e_L = prob_e_s[(e_L, s_H)] + prob_e_s[(e_L, s_L)]

    # Conditional probability P(s|e) = P(e,s) / P(e)
    P_sH_given_eH = prob_e_s[(e_H, s_H)] / P_e_H
    P_sL_given_eH = prob_e_s[(e_H, s_L)] / P_e_H

    # Conditional expectation E[s|e] = s_H * P(s_H|e) + s_L * P(s_L|e)
    # The high effort e_H is chosen when the worker gets a signal for the high state.
    E_s_given_eH = s_H * P_sH_given_eH + s_L * P_sL_given_eH

    print("The worker's optimal effort `e` is determined by the equation: e = beta * E[s | signal].")
    print(f"The conditional expectation of the state, given the high effort signal, is calculated as:")
    print(f"E[s | e={e_H}] = s_H * P(s={s_H}|e={e_H}) + s_L * P(s={s_L}|e={e_H})")
    print(f"E[s | e={e_H}] = {s_H} * {P_sH_given_eH:.4f} + {s_L} * {P_sL_given_eH:.4f} = {E_s_given_eH}\n")

    # Step 4: Determine beta from the worker's incentive compatibility constraint
    # e_H = beta * E[s | e_H]
    beta = e_H / E_s_given_eH
    print("We can now solve for the contract parameter beta using the high effort case:")
    print(f"Equation: {e_H} = beta * {E_s_given_eH}")
    print(f"beta = {e_H} / {E_s_given_eH}")
    print(f"beta = {beta}\n")
    
    # Step 5 & 6: Relate beta to p and find the final answer
    # For a risk-neutral agent, the firm's profit is maximized by setting beta = p.
    p = beta
    print("In this model with a risk-neutral employee, the firm maximizes profit by setting beta equal to the output price p.")
    print("Final Equation: p = beta")
    print(f"Therefore, the value of p is {p}.")

    return p

if __name__ == '__main__':
    final_p = solve_for_p()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{final_p}>>>")