import collections

def solve_price():
    """
    Calculates the price 'p' based on the given probability distribution
    of effort and output.
    """
    # Step 1: Define the given data and deduce states and efforts.
    # P(e, y) is given
    prob_e_y = {
        (22, 132): 0.4375,
        (22, 44): 0.0625,
        (10, 60): 0.0625,
        (10, 20): 0.4375
    }

    print("Step 1: Analyzing the given data P(e, y)")
    print("=========================================")
    print("Given probability distribution P(e, y):")
    for (e, y), prob in prob_e_y.items():
        print(f"P(e={e}, y={y}) = {prob}")

    # Deduce states s = y/e
    states = sorted(list(set([y / e for (e, y) in prob_e_y.keys()])))
    s_L, s_H = states[0], states[1]
    efforts = sorted(list(set([e for (e, y) in prob_e_y.keys()])), reverse=True)
    e_H, e_L = efforts[0], efforts[1]
    
    print(f"\nDeduced states of the world: s_L = {s_L}, s_H = {s_H}")
    print(f"Effort levels: e_H = {e_H} (high), e_L = {e_L} (low)")

    # Create the joint probability distribution P(e, s)
    prob_e_s = collections.defaultdict(float)
    for (e, y), prob in prob_e_y.items():
        s = y / e
        prob_e_s[(e, s)] = prob
    
    print("\nConverted to joint probability distribution P(e, s):")
    for (e, s), prob in prob_e_s.items():
        print(f"P(e={e}, s={s}) = {prob}")

    # Step 2: Calculate marginal and conditional probabilities.
    print("\nStep 2: Calculating conditional probabilities P(s|e)")
    print("===================================================")
    
    # Marginal probabilities of effort P(e)
    prob_e = collections.defaultdict(float)
    for (e, s), prob in prob_e_s.items():
        prob_e[e] += prob

    print(f"P(e={e_H}) = {prob_e_s[(e_H, s_H)]} + {prob_e_s[(e_H, s_L)]} = {prob_e[e_H]}")
    print(f"P(e={e_L}) = {prob_e_s[(e_L, s_H)]} + {prob_e_s[(e_L, s_L)]} = {prob_e[e_L]}")
    
    # Conditional probabilities P(s|e) = P(e,s) / P(e)
    prob_s_given_e = collections.defaultdict(float)
    for (e, s), prob in prob_e_s.items():
        prob_s_given_e[(s, e)] = prob / prob_e[e]

    print("\nConditional probabilities P(s|e):")
    print(f"P(s={s_H}|e={e_H}) = {prob_e_s[(e_H, s_H)]} / {prob_e[e_H]} = {prob_s_given_e[(s_H, e_H)]}")
    print(f"P(s={s_L}|e={e_H}) = {prob_e_s[(e_H, s_L)]} / {prob_e[e_H]} = {prob_s_given_e[(s_L, e_H)]}")
    print(f"P(s={s_H}|e={e_L}) = {prob_e_s[(e_L, s_H)]} / {prob_e[e_L]} = {prob_s_given_e[(s_H, e_L)]}")
    print(f"P(s={s_L}|e={e_L}) = {prob_e_s[(e_L, s_L)]} / {prob_e[e_L]} = {prob_s_given_e[(s_L, e_L)]}")

    # Step 3: Calculate conditional expectations E[s|e].
    print("\nStep 3: Calculating conditional expectations E[s|e]")
    print("==================================================")
    
    E_s_given_eH = s_H * prob_s_given_e[(s_H, e_H)] + s_L * prob_s_given_e[(s_L, e_H)]
    E_s_given_eL = s_H * prob_s_given_e[(s_H, e_L)] + s_L * prob_s_given_e[(s_L, e_L)]
    
    print(f"E[s|e={e_H}] = {s_H}*P(s={s_H}|e={e_H}) + {s_L}*P(s={s_L}|e={e_H})")
    print(f"E[s|e={e_H}] = {s_H}*{prob_s_given_e[(s_H, e_H)]} + {s_L}*{prob_s_given_e[(s_L, e_H)]} = {E_s_given_eH}")
    
    print(f"E[s|e={e_L}] = {s_H}*P(s={s_H}|e={e_L}) + {s_L}*P(s={s_L}|e={e_L})")
    print(f"E[s|e={e_L}] = {s_H}*{prob_s_given_e[(s_H, e_L)]} + {s_L}*{prob_s_given_e[(s_L, e_L)]} = {E_s_given_eL}")

    # Step 4: Determine beta from the employee's optimal effort rule e = beta * E[s|e].
    print("\nStep 4: Calculating beta (β)")
    print("============================")
    
    beta = e_H / E_s_given_eH
    print(f"From the high effort level, e_H = β * E[s|e_H]")
    print(f"{e_H} = β * {E_s_given_eH}")
    print(f"β = {e_H} / {E_s_given_eH} = {beta}")

    # Step 5: State the theoretical result that p = beta.
    print("\nStep 5: Final Conclusion")
    print("=========================")
    print("In this principal-agent model, the firm's profit-maximizing choice for the contract parameter β is equal to the output price p.")
    print("Therefore, p = β.")
    
    # Step 6: Final Answer
    p = beta
    print(f"The value of p is {p}.")
    
    return p

if __name__ == '__main__':
    price = solve_price()
    print(f"\n<<<__{price}__>>>")
