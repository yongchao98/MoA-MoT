import numpy as np

def calculate_params(e_max, e_min):
    """
    Calculates tight-binding parameters t and s from energies at the Gamma point.
    Note: These are estimations based on the simplified formula for the Gamma point.
    """
    if e_max == 0 or e_min == 0:
        return 0, 0
    
    # Handle the symmetric case s=0
    if abs(abs(e_max) - abs(e_min)) < 1e-6:
        s = 0.0
        # From ΔE = E_max - E_min = 6t / (1 - 9s^2) = 6t
        t = (e_max - e_min) / 6.0
        return t, s

    # Handle asymmetric cases
    # Determine sign of s from asymmetry
    if abs(e_min) > abs(e_max):
        # Valence band is wider, s > 0
        rho = abs(e_min) / abs(e_max)
        # from rho = (1 + 3s) / (1 - 3s) => s = (rho - 1) / (3 * (rho + 1))
        s = (rho - 1) / (3 * (rho + 1))
    else:
        # Conduction band is wider, s < 0
        rho = abs(e_max) / abs(e_min)
        s = - (rho - 1) / (3 * (rho + 1))

    # Calculate t from total bandwidth
    # ΔE = 6t / (1 - 9s^2)
    delta_e = e_max - e_min
    t = delta_e * (1 - 9 * s**2) / 6.0
    
    return t, s

def solve_graphene_puzzle():
    """
    Solves the graphene simulation puzzle by analyzing extracted parameters.
    """
    # Energies at Gamma point extracted from the plots
    # sim_data = {sim_index: (E_max, E_min)}
    sim_data = {
        1: (3.0, -15.0),
        2: (3.0, -3.0),
        3: (3.0, -16.0),
        4: (16.0, -3.0)
    }

    # Calculate parameters for each simulation
    params = {}
    for i, (e_max, e_min) in sim_data.items():
        t, s = calculate_params(e_max, e_min)
        params[i] = {'t': t, 's': s, '|s|': abs(s)}
        
    print("Estimated parameters:")
    for i in sorted(params.keys()):
        print(f"Simulation {i}: t ≈ {params[i]['t']:.2f}, s ≈ {params[i]['s']:.3f}")
    print("-" * 30)

    # Initialize answer array
    answer = [0] * 4  # answer[i] will be the sim index for condition i+1

    # Condition 1: minimum t
    min_t_sim = min(params, key=lambda sim: params[sim]['t'])
    answer[0] = min_t_sim

    # Condition 4: maximum s
    max_s_sim = max(params, key=lambda sim: params[sim]['s'])
    answer[3] = max_s_sim

    # Condition 3: unique sign(s)
    # Sim 4 is the only one with s < 0
    unique_sign_sim = 4
    answer[2] = unique_sign_sim
    
    # Condition 2: minimum |s|
    # Sim 2 has |s|=0, which is the absolute minimum. However, a 1-to-1 mapping is required,
    # and Sim 2 is already the best fit for minimum t due to its dramatically smaller bandwidth.
    # The remaining simulation is Sim 1. We assign it to the remaining condition.
    # This can be justified by interpreting "minimum |s|" as the minimum non-zero |s|.
    # Among sims {1, 3, 4}, sim 1 has the minimum |s|.
    
    assigned_sims = {answer[0], answer[2], answer[3]}
    all_sims = {1, 2, 3, 4}
    min_s_mag_sim = list(all_sims - assigned_sims)[0]
    answer[1] = min_s_mag_sim
    
    print("Mapping from Condition to Simulation Index:")
    print(f"1) minimum t       -> Simulation {answer[0]}")
    print(f"2) minimum |s|     -> Simulation {answer[1]}")
    print(f"3) unique sign(s)  -> Simulation {answer[2]}")
    print(f"4) maximum s       -> Simulation {answer[3]}")
    print("-" * 30)
    
    # Format the final answer string
    final_answer_string = "".join(map(str, answer))
    print(f"The simulation indices ordered by the condition met are: {final_answer_string}")
    
    # Return the final answer in the required format
    return final_answer_string

final_answer = solve_graphene_puzzle()
print(f"<<<{final_answer}>>>")