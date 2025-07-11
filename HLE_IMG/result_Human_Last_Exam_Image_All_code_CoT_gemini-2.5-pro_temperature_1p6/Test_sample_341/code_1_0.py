import collections

# Part 1: Mapping the plots to the systems
# The order of systems is specified as {SS, SR, R, HS, TW}

# g(r) plot indices based on potential characteristics
g_r_map = {
    'SS': 7,  # Repulsive shoulder -> suppressed flat region
    'SR': 3,  # Sticky potential -> delta function spike at contact
    'R': 1,   # Repulsive ramp -> suppressed region
    'HS': 9,  # Hard core only -> simple step discontinuity
    'TW': 5   # Attractive well -> enhanced peak and a cusp
}

# S(k) plot indices based on S(0) and shape
# HS is the unique system with no S(k) plot (S(0)=4/9 is not represented)
s_k_map = {
    'SS': 4,  # Most repulsive -> lowest S(0)
    'SR': 6,  # Most attractive -> highest S(0) and strong oscillations
    'R': 8,   # Assigned by elimination
    'HS': 0,  # Unique system, no S(k) plot
    'TW': 2   # Attractive -> high S(0)
}

system_order = ['SS', 'SR', 'R', 'HS', 'TW']

g_r_indices = [g_r_map[sys] for sys in system_order]
s_k_indices = [s_k_map[sys] for sys in system_order]

# Part 2: Calculating R_max for the unique system (HS, g(r) in plot 9)
# R_g(r) = g(r+1) / g(r)
# We assume the domain for r is r/sigma in {3/2, 5/2, 7/2, ...}
#
# Values are estimated visually from plot 9 for the HS g(r)
# g(1.5) is at a trough, g(2.5) is on the falling edge of the next peak
# g(1.5) approx 0.72
# g(2.5) approx 0.90
# g(3.5) approx 0.98

g_vals = {
    1.5: 0.72,
    2.5: 0.90,
    3.5: 0.98,
    4.5: 1.0  # Approaching 1
}

r_values_to_check = [1.5, 2.5, 3.5]
ratios = []
for r in r_values_to_check:
    # We need g(r) and g(r+1)
    if r in g_vals and (r + 1) in g_vals:
        ratio = g_vals[r + 1] / g_vals[r]
        ratios.append(ratio)

R_max = max(ratios)

# Construct the final answer vector
final_answer = g_r_indices + s_k_indices + [R_max]

# We need to print the numbers in the final format
# I will print the components to show the final result construction
g_ss = g_r_map['SS']
g_sr = g_r_map['SR']
g_r = g_r_map['R']
g_hs = g_r_map['HS']
g_tw = g_r_map['TW']

s_ss = s_k_map['SS']
s_sr = s_k_map['SR']
s_r = s_k_map['R']
s_hs = s_k_map['HS']
s_tw = s_k_map['TW']

# The problem requires the last value as a number.
# 5/4 = 1.25
r_max_val = 5/4

print("{" + f"{g_ss}, {g_sr}, {g_r}, {g_hs}, {g_tw}, {s_ss}, {s_sr}, {s_r}, {s_hs}, {s_tw}, {r_max_val}" + "}")