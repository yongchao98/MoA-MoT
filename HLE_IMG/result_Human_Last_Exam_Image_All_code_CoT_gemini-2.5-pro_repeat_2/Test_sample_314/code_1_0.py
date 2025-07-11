import numpy as np

def calculate_params_from_gamma(E_plus_gamma, E_minus_gamma):
    """
    Calculates t and s from energy values at the Gamma point.
    E_plus = 3t / (1 + 3s)
    E_minus = -3t / (1 - 3s)
    """
    if E_plus_gamma == -E_minus_gamma: # Symmetric case
        s = 0
        if E_plus_gamma > 0:
            t = E_plus_gamma / 3.0
        else: # Should not happen for E_plus
            t = 0
        return t, s

    # Asymmetric case
    # E_plus / E_minus = -(1 - 3s) / (1 + 3s)
    ratio = -E_plus_gamma / E_minus_gamma
    # ratio * (1 + 3s) = 1 - 3s
    # ratio + 3s*ratio = 1 - 3s
    # 3s(ratio + 1) = 1 - ratio
    s = (1 - ratio) / (3 * (ratio + 1))
    
    # 3t = E_plus * (1 + 3s)
    t = E_plus_gamma * (1 + 3 * s) / 3.0
    return t, s

def calculate_params_from_gamma_M(E_minus_gamma, E_plus_M):
    """
    Calculates t and s from E-(Gamma) and E+(M).
    -E_minus_gamma = 3t / (1 - 3s)
    E_plus_M = t / (1 + s)
    From the second eq: t = E_plus_M * (1 + s)
    Substitute into first: -E_minus_gamma = 3 * E_plus_M * (1+s) / (1-3s)
    -E_minus_gamma * (1-3s) = 3 * E_plus_M * (1+s)
    -E_minus_gamma + 3s*E_minus_gamma = 3*E_plus_M + 3s*E_plus_M
    s * (3*E_minus_gamma - 3*E_plus_M) = 3*E_plus_M + E_minus_gamma
    s = (3*E_plus_M + E_minus_gamma) / (3*E_minus_gamma - 3*E_plus_M)
    """
    s = (3 * E_plus_M + E_minus_gamma) / (3 * E_minus_gamma - 3 * E_plus_M)
    t = E_plus_M * (1 + s)
    return t, s

# --- Parameter Extraction from Plots ---

# Simulation 1: Read from 2D plot K-Gamma-M-K'
# E-(Gamma) = -15 eV, E+(M) = 2.5 eV
t1, s1 = calculate_params_from_gamma_M(-15, 2.5)

# Simulation 2: Symmetric 3D plot
# E_max = 10, E_min = -10. This corresponds to E_plus(Gamma) and E_minus(Gamma) for s=0.
t2, s2 = calculate_params_from_gamma(10, -10)

# Simulation 3: Asymmetric 3D plot, bottom-heavy
# E_plus(Gamma) ~ 5, E_minus(Gamma) ~ -15
t3, s3 = calculate_params_from_gamma(5, -15)

# Simulation 4: Asymmetric 3D plot, top-heavy
# E_plus(Gamma) ~ 15, E_minus(Gamma) ~ -5
t4, s4 = calculate_params_from_gamma(15, -5)

simulations = {
    1: {'t': t1, 's': s1},
    2: {'t': t2, 's': s2},
    3: {'t': t3, 's': s3},
    4: {'t': t4, 's': s4},
}

print("Calculated Tight-Binding Parameters:")
for i in sorted(simulations.keys()):
    print(f"Simulation {i}: t = {simulations[i]['t']:.3f}, s = {simulations[i]['s']:.3f}")

# --- Logic to assign simulations to conditions ---

conditions = {
    "min_t": None,
    "min_abs_s": None,
    "unique_sign_s": None,
    "max_s": None,
}

# Create a list of available simulations
available_sims = list(simulations.keys())
result = [0, 0, 0, 0] # Index will be condition number - 1

# Condition 2: minimum |s|
min_abs_s_val = min(abs(p['s']) for p in simulations.values())
min_s_sims = [i for i, p in simulations.items() if abs(p['s']) == min_abs_s_val]
if len(min_s_sims) == 1:
    sim_idx = min_s_sims[0]
    conditions["min_abs_s"] = sim_idx
    result[1] = sim_idx
    available_sims.remove(sim_idx)

# Condition 4: maximum s
max_s_val = max(p['s'] for p in simulations.values())
max_s_sims = [i for i, p in simulations.items() if p['s'] == max_s_val]
if len(max_s_sims) == 1:
    sim_idx = max_s_sims[0]
    conditions["max_s"] = sim_idx
    if sim_idx in available_sims:
        result[3] = sim_idx
        available_sims.remove(sim_idx)

# Condition 3: unique sign(s)
signs = [np.sign(p['s']) for p in simulations.values()]
sign_counts = {sign: signs.count(sign) for sign in set(signs)}
unique_signs = [sign for sign, count in sign_counts.items() if count == 1]
# We are looking for the unique non-zero sign, which is -1.
if -1 in unique_signs:
    unique_sign_sims = [i for i, p in simulations.items() if np.sign(p['s']) == -1]
    if len(unique_sign_sims) == 1:
        sim_idx = unique_sign_sims[0]
        conditions["unique_sign_s"] = sim_idx
        if sim_idx in available_sims:
            result[2] = sim_idx
            available_sims.remove(sim_idx)

# Condition 1: minimum t is assigned to the last remaining simulation
if len(available_sims) == 1:
    sim_idx = available_sims[0]
    conditions["min_t"] = sim_idx
    result[0] = sim_idx
    available_sims.remove(sim_idx)

print("\nMapping of conditions to simulations:")
print(f"1) minimum t:       Simulation {result[0]}")
print(f"2) minimum |s|:     Simulation {result[1]}")
print(f"3) unique sign(s):  Simulation {result[2]}")
print(f"4) maximum s:       Simulation {result[3]}")

final_answer_string = "".join(map(str, result))
print(f"\nThe final ordered answer is the sequence of simulation indices: {final_answer_string}")

# Final Answer Block
print("\nFinal Answer in required format:")
print(f"<<<{final_answer_string}>>>")