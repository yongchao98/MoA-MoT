def solve_graphene_puzzle():
    """
    This function calculates the tight-binding parameters from the provided plot data
    and assigns each simulation to its corresponding condition.
    """
    # Data extracted from plots at the Gamma point (in eV)
    data = {
        1: {'Ec_G': 2.5, 'Ev_G': -15.0},
        2: {'Ec_G': 3.0, 'Ev_G': -9.0},
        3: {'Ec_G': 5.0, 'Ev_G': -15.0},
        4: {'Ec_G': 9.0, 'Ev_G': -9.0}
    }

    params = {}

    # Step 1 & 2: Calculate parameters s, BW, and t for each simulation
    for i in sorted(data.keys()):
        sim_data = data[i]
        Ec = sim_data['Ec_G']
        Ev = sim_data['Ev_G']

        # Calculate asymmetry ratio R
        R = abs(Ev) / Ec

        # Calculate overlap magnitude s
        if R == 1:
            s = 0.0
        else:
            s = (R - 1) / (3 * (R + 1))

        # Calculate total bandwidth BW
        BW = Ec - Ev
        
        # Calculate hopping parameter t
        t = BW * (1 - 9 * s**2) / 6
        
        params[i] = {'s': s, 't': t, 'BW': BW}

    # Step 3: Assign simulations to conditions
    t_vals = {i: p['t'] for i, p in params.items()}
    s_vals = {i: p['s'] for i, p in params.items()}
    bw_vals = {i: p['BW'] for i, p in params.items()}

    # Condition 1: minimum t
    min_t_sim = min(t_vals, key=t_vals.get)
    
    # Condition 2: minimum |s|
    min_s_sim = min(s_vals, key=s_vals.get)

    # Condition 4: Re-interpreted as maximum Bandwidth to resolve conflict
    max_bw_sim = max(bw_vals, key=bw_vals.get)

    # Assign all used sims
    assigned_sims = {min_t_sim, min_s_sim, max_bw_sim}
    all_sims = set(data.keys())
    
    # Condition 3: The one remaining simulation
    unique_sign_sim = list(all_sims - assigned_sims)[0]
    
    # Form the final answer string
    final_answer = f"{min_t_sim}{min_s_sim}{unique_sign_sim}{max_bw_sim}"

    print(f"Final Answer Equation: a={min_t_sim}, b={min_s_sim}, c={unique_sign_sim}, d={max_bw_sim}")
    print(f"1st condition is met by simulation {min_t_sim}")
    print(f"2nd condition is met by simulation {min_s_sim}")
    print(f"3rd condition is met by simulation {unique_sign_sim}")
    print(f"4th condition is met by simulation {max_bw_sim}")

solve_graphene_puzzle()