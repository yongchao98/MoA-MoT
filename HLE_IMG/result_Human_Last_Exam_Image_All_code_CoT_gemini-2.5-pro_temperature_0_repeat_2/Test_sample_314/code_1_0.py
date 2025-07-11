import numpy as np

def calculate_params(Ec, Ev):
    """
    Calculates tight-binding parameters t and s from conduction and valence band energies at the Gamma point.
    """
    # Determine sign of s from asymmetry
    if abs(Ev) > abs(Ec):
        # s > 0
        # Ec/(-Ev) = (1-3s)/(1+3s) => s = (-Ev-Ec)/(3*(Ec-Ev))
        s = (-Ev - Ec) / (3 * (Ec - Ev))
        # t = Ec*(1+3s)/3
        t = Ec * (1 + 3 * s) / 3
    else:
        # s < 0, let s = -|s|
        # Ec/(-Ev) = (1+3|s|)/(1-3|s|) => |s| = (-Ev-Ec)/(3*(Ev-Ec))
        s_mag = (-Ev - Ec) / (3 * (Ev - Ec))
        s = -s_mag
        # t = Ec*(1-3|s|)/3
        t = Ec * (1 - 3 * s_mag) / 3
    return t, s

def solve_puzzle():
    """
    Solves the puzzle by calculating parameters and matching them to conditions.
    """
    # Step 2: Extract data from plots
    energies = {
        1: {'Ec': 3, 'Ev': -15},
        2: {'Ec': 3, 'Ev': -10},
        3: {'Ec': 5, 'Ev': -15},
        4: {'Ec': 15, 'Ev': -10},
    }

    # Step 3: Calculate t and s for each simulation
    params = {}
    for i in range(1, 5):
        t, s = calculate_params(energies[i]['Ec'], energies[i]['Ev'])
        params[i] = {'t': t, 's': s}
        # print(f"Simulation {i}: t = {t:.3f}, s = {s:.3f}")

    # Step 4: Apply conditions to find the mapping
    
    # Condition 3: unique sign(s)
    # Find the simulation with a unique sign for s
    s_signs = [np.sign(params[i]['s']) for i in range(1, 5)]
    unique_sign_sim = 0
    for i in range(1, 5):
        if s_signs.count(np.sign(params[i]['s'])) == 1:
            unique_sign_sim = i
            break
    
    # The remaining simulations and their indices
    remaining_indices = [i for i in range(1, 5) if i != unique_sign_sim]
    
    # Apply other conditions to the remaining simulations
    
    # Condition 1: minimum t
    min_t_val = float('inf')
    min_t_sim = 0
    for i in remaining_indices:
        if params[i]['t'] < min_t_val:
            min_t_val = params[i]['t']
            min_t_sim = i

    # Condition 2: minimum |s|
    min_s_val = float('inf')
    min_s_sim = 0
    for i in remaining_indices:
        if abs(params[i]['s']) < min_s_val:
            min_s_val = abs(params[i]['s'])
            min_s_sim = i
            
    # Condition 4: maximum s
    max_s_val = float('-inf')
    max_s_sim = 0
    for i in remaining_indices:
        if params[i]['s'] > max_s_val:
            max_s_val = params[i]['s']
            max_s_sim = i

    # Step 5: Assemble the final ordered list
    # The order is based on the condition number (1, 2, 3, 4)
    result = [min_t_sim, min_s_sim, unique_sign_sim, max_s_sim]
    
    print("Based on the analysis of the band structures:")
    print(f"1. The simulation with the minimum hopping parameter (t) is: {min_t_sim}")
    print(f"2. The simulation with the minimum overlap magnitude (|s|) among those with the majority sign is: {min_s_sim}")
    print(f"3. The simulation with the unique overlap sign (sign(s)) is: {unique_sign_sim}")
    print(f"4. The simulation with the maximum overlap (s) among those with the majority sign is: {max_s_sim}")
    
    print("\nThe final answer, with simulation indices ordered by the condition met, is:")
    # Print the result as a single number string
    final_answer = "".join(map(str, result))
    print(final_answer)

solve_puzzle()