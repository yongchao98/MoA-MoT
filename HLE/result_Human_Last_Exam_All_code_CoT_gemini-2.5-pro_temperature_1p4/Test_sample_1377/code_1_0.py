import math

def calculate_rates():
    """
    Calculates the optimal growth rate W_star and the actual growth rate W.
    """
    # True probabilities
    p = {'A': 0.5, 'B': 0.25, 'C': 0.25}
    # Mistaken probabilities
    q = {'A': 0.25, 'B': 0.5, 'C': 0.25}
    # Net odds (b-to-1)
    b = {'A': 4, 'B': 3, 'C': 3}

    print("Step 1: Calculate the optimal strategy and growth rate (W*)")

    # For the rational bettor with true probabilities p
    # Calculate single-bet Kelly fractions
    f = {}
    for i in p:
        if p[i] * b[i] - (1 - p[i]) > 0:
            f[i] = p[i] - (1 - p[i]) / b[i]
        else:
            f[i] = 0

    print(f"Optimal Kelly fractions for single bets: f_A={f['A']:.4f}, f_B={f['B']:.4f}, f_C={f['C']:.4f}")

    # Calculate growth rate for each potential single bet
    W_star_options = {}
    if f['A'] > 0:
        W_star_options['A'] = p['A'] * math.log(1 + f['A'] * b['A']) + (1 - p['A']) * math.log(1 - f['A'])
    else:
        W_star_options['A'] = 0

    if f['B'] > 0:
        W_star_options['B'] = p['B'] * math.log(1 + f['B'] * b['B']) + (1 - p['B']) * math.log(1 - f['B'])
    else:
        W_star_options['B'] = 0
    
    if f['C'] > 0:
        W_star_options['C'] = p['C'] * math.log(1 + f['C'] * b['C']) + (1 - p['C']) * math.log(1 - f['C'])
    else:
        W_star_options['C'] = 0

    # The optimal strategy is to pick the best of these single bets
    W_star = max(W_star_options.values())
    rational_choice = max(W_star_options, key=W_star_options.get)
    rational_fraction = f[rational_choice]

    print(f"Growth rate if betting on A: {W_star_options['A']:.4f}")
    print(f"Growth rate if betting on B: {W_star_options['B']:.4f}")
    print(f"Growth rate if betting on C: {W_star_options['C']:.4f}")
    print(f"The rational bettor chooses to bet on {rational_choice} with fraction {rational_fraction:.4f}.")
    print(f"W* = {W_star:.6f}\n")


    print("Step 2: Calculate the mistaken bettor's strategy")
    
    # For the mistaken bettor with probabilities q
    g = {}
    for i in q:
        if q[i] * b[i] - (1 - q[i]) > 0:
            g[i] = q[i] - (1 - q[i]) / b[i]
        else:
            g[i] = 0
            
    print(f"Mistaken Kelly fractions for single bets: g_A={g['A']:.4f}, g_B={g['B']:.4f}, g_C={g['C']:.4f}")
    
    # Calculate perceived growth rate for each potential bet
    W_prime_options = {}
    if g['A'] > 0:
        W_prime_options['A'] = q['A'] * math.log(1 + g['A'] * b['A']) + (1 - q['A']) * math.log(1 - g['A'])
    else:
        W_prime_options['A'] = 0

    if g['B'] > 0:
        W_prime_options['B'] = q['B'] * math.log(1 + g['B'] * b['B']) + (1 - q['B']) * math.log(1 - g['B'])
    else:
        W_prime_options['B'] = 0
        
    if g['C'] > 0:
        W_prime_options['C'] = q['C'] * math.log(1 + g['C'] * b['C']) + (1 - q['C']) * math.log(1 - g['C'])
    else:
        W_prime_options['C'] = 0

    mistaken_choice = max(W_prime_options, key=W_prime_options.get)
    mistaken_fraction = g[mistaken_choice]
    
    print(f"Perceived growth if betting on A: {W_prime_options['A']:.4f}")
    print(f"Perceived growth if betting on B: {W_prime_options['B']:.4f}")
    print(f"Perceived growth if betting on C: {W_prime_options['C']:.4f}")
    print(f"The mistaken bettor chooses to bet on {mistaken_choice} with fraction {mistaken_fraction:.4f}.\n")


    print("Step 3: Calculate the actual growth rate (W) for the mistaken bettor")

    # The actual growth is based on the mistaken strategy but with true probabilities
    # The mistaken bettor bets fraction g_B on horse B.
    # We calculate the growth rate using the true probabilities p.
    if mistaken_choice == 'A':
        wealth_if_A_wins = 1 + mistaken_fraction * b['A']
        wealth_if_B_wins = 1 - mistaken_fraction
        wealth_if_C_wins = 1 - mistaken_fraction
    elif mistaken_choice == 'B':
        wealth_if_A_wins = 1 - mistaken_fraction
        wealth_if_B_wins = 1 + mistaken_fraction * b['B']
        wealth_if_C_wins = 1 - mistaken_fraction
    elif mistaken_choice == 'C':
        wealth_if_A_wins = 1 - mistaken_fraction
        wealth_if_B_wins = 1 - mistaken_fraction
        wealth_if_C_wins = 1 + mistaken_fraction * b['C']

    W = (p['A'] * math.log(wealth_if_A_wins) +
         p['B'] * math.log(wealth_if_B_wins) +
         p['C'] * math.log(wealth_if_C_wins))

    print(f"Actual growth rate for mistaken strategy (W) = {p['A']:.2f}*log({wealth_if_A_wins:.4f}) + "
          f"{p['B']:.2f}*log({wealth_if_B_wins:.4f}) + {p['C']:.2f}*log({wealth_if_C_wins:.4f})")
    print(f"W = {W:.6f}\n")
    
    
    print("Step 4: Calculate W* - W")
    
    # W_star = 0.5 * math.log(1+0.375*4) + 0.5 * math.log(1-0.375)
    # W = 0.5 * math.log(1-1/3) + 0.25 * math.log(1+(1/3)*3) + 0.25 * math.log(1-1/3)
    diff = W_star - W
    
    print(f"The final equation is (W*) - (W):")
    print(f"({p['A']}*log(1 + {f['A']}*{b['A']}) + {1-p['A']}*log(1 - {f['A']})) - ")
    print(f"({p['A']}*log(1-{g['B']:.4f}) + {p['B']}*log(1 + {g['B']:.4f}*{b['B']}) + {p['C']}*log(1-{g['B']:.4f}))")
    print(f"= ({W_star:.6f}) - ({W:.6f})")
    print(f"= {diff:.6f}")
    
    return diff

# Run the calculation and print the final answer
result = calculate_rates()
print(f"\nFinal Answer: {result:.6f}")