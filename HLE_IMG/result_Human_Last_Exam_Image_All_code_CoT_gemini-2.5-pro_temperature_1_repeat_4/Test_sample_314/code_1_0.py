import numpy as np

def calculate_parameters(E_max, E_min, s_is_positive):
    """
    Calculates tight-binding parameters t and s from energy values at the Gamma point.
    f_gamma is the geometric factor at the Gamma point, which is 3.
    """
    f_gamma = 3.0
    
    if s_is_positive:
        # E_upper = t*f / (1 - s*f)  =>  t = E_upper * (1 - s*f) / f
        # E_lower = -t*f / (1 + s*f) =>  t = -E_lower * (1 + s*f) / f
        # E_upper * (1 - s*f) = -E_lower * (1 + s*f)
        # E_upper - s*f*E_upper = -E_lower - s*f*E_lower
        # E_upper + E_lower = s*f*(E_upper - E_lower)
        # s = (E_upper + E_lower) / (f * (E_upper - E_lower))
        s = (E_max + E_min) / (f_gamma * (E_max - E_min))
        t = E_max * (1 - s * f_gamma) / f_gamma
    else: # s is negative
        # E_upper = t*f / (1 + |s|*f) => t = E_upper * (1 + |s|*f) / f
        # E_lower = -t*f / (1 - |s|*f) => t = -E_lower * (1 - |s|*f) / f
        # E_upper * (1 + |s|*f) = -E_lower * (1 - |s|*f)
        # E_upper + |s|*f*E_upper = -E_lower + |s|*f*E_lower
        # E_upper + E_lower = |s|*f*(E_lower - E_upper)
        # |s| = (E_upper + E_lower) / (f*(E_lower - E_upper))
        s_mag = (E_max + E_min) / (f_gamma * (E_min - E_max))
        s = -s_mag
        t = E_max * (1 + s_mag * f_gamma) / f_gamma
        
    return t, s

# Data extracted from plots
simulations = {
    1: {'E_max': 2.5, 'E_min': -15, 's_is_positive': False},
    2: {'E_max': 3.0, 'E_min': -10, 's_is_positive': False},
    3: {'E_max': 5.0, 'E_min': -15, 's_is_positive': False},
    4: {'E_max': 15.0, 'E_min': -5, 's_is_positive': True}
}

# Calculate parameters for each simulation
params = {}
for i, data in simulations.items():
    t, s = calculate_parameters(data['E_max'], data['E_min'], data['s_is_positive'])
    params[i] = {'t': t, 's': s, '|s|': abs(s)}
    print(f"Simulation {i}: t = {t:.3f}, s = {s:.3f}, |s| = {abs(s):.3f}")

# Determine which simulation meets each condition
t_values = {i: p['t'] for i, p in params.items()}
s_values = {i: p['s'] for i, p in params.items()}
s_mag_values = {i: p['|s|'] for i, p in params.items()}

# Initialize answer array and set of used indices
ans = [0, 0, 0, 0]
used_indices = set()

# Condition 3: unique sign(s)
# This is the most distinct qualitative feature. Sim 4 is the only one with s>0.
ans[2] = 4
used_indices.add(4)
print(f"\nCondition 3 (unique sign(s)) is met by simulation {ans[2]}.")

# Condition 1: minimum t
# Find the simulation with the minimum t value among the remaining ones.
min_t = float('inf')
min_t_idx = -1
for i in range(1, 5):
    if i not in used_indices:
        if t_values[i] < min_t:
            min_t = t_values[i]
            min_t_idx = i
ans[0] = min_t_idx
used_indices.add(min_t_idx)
print(f"Condition 1 (minimum t) is met by simulation {ans[0]}.")

# Condition 2: minimum |s|
# Find the simulation with the minimum |s| value among the remaining ones.
min_s_mag = float('inf')
min_s_mag_idx = -1
for i in range(1, 5):
    if i not in used_indices:
        # Use a small tolerance for comparison due to reading from plots
        if s_mag_values[i] < min_s_mag - 1e-9:
            min_s_mag = s_mag_values[i]
            min_s_mag_idx = i
ans[1] = min_s_mag_idx
used_indices.add(min_s_mag_idx)
print(f"Condition 2 (minimum |s|) is met by simulation {ans[1]}.")

# Condition 4: maximum s
# The last remaining index is assigned to the last condition to form the required permutation.
# There is a conflict in the problem statement as 'max s' also points to Sim 4.
# We resolve this by assigning the last available index.
for i in range(1, 5):
    if i not in used_indices:
        ans[3] = i
        used_indices.add(i)
print(f"Condition 4 (maximum s) is assigned to the last available simulation, {ans[3]}, to form a permutation.")


final_answer_str = "".join(map(str, ans))
print(f"\nThe final ordered list of simulation indices is: {ans[0]}, {ans[1]}, {ans[2]}, {ans[3]}")
print(f"Final answer string: {final_answer_str}")
# The final answer is the concatenation of these numbers.
# We must print the equation in the final result.
# I will print the final answer in the requested format.
print(f"Final Answer: {ans[0]}{ans[1]}{ans[2]}{ans[3]}")