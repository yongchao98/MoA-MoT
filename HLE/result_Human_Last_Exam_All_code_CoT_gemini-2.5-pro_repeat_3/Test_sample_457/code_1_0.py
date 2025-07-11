import numpy as np

# Step 1: Define the provided data
# Data as two NumPy arrays: relative pressure (P/P0) and adsorbed volume (V)
p_p0_raw = np.array([
    0.000097778, 0.000490286, 0.001028258, 0.005256671, 0.010218832, 0.053231858,
    0.098896253, 0.14976513, 0.200004704, 0.233517997, 0.266801197, 0.301984569,
    0.325449461, 0.35054002, 0.376514193, 0.400208966, 0.426400846, 0.450079581,
    0.476335676, 0.499897655, 0.525467416, 0.551884881, 0.575255893, 0.600463037,
    0.625966848, 0.651702807, 0.676434592, 0.701749025, 0.726670645, 0.751424402,
    0.776687236, 0.801317436, 0.826598164, 0.85144159, 0.876516006, 0.901802511,
    0.926576152, 0.96028434, 0.979849853, 0.981764776, 0.989523087, 0.999369551,
    0.984607125, 0.967573883, 0.944626529, 0.927201724, 0.90187396, 0.877154076,
    0.851603484, 0.826517066, 0.802151817, 0.777231389, 0.751936033, 0.727508526,
    0.702376073, 0.676883722, 0.650369795, 0.62669699, 0.601288503, 0.576203003,
    0.550292959, 0.525392205, 0.501925833, 0.475695455, 0.447680802, 0.423654664,
    0.382079135, 0.356158655, 0.331104613, 0.325044783, 0.299552276
])

v_ads_raw = np.array([
    34.912, 52.8679, 63.2559, 90.9391, 104.0776, 144.7389, 165.4765, 182.3897,
    196.5799, 205.2201, 213.3343, 221.7435, 227.2067, 233.1827, 239.264, 244.933,
    251.3238, 257.1772, 263.9726, 270.2935, 277.4911, 285.6177, 293.3731, 302.0746,
    312.2637, 323.6833, 335.6003, 348.5889, 362.0965, 376.0429, 390.2814, 404.5596,
    417.7411, 430.0732, 440.2639, 446.8344, 450.299, 454.6433, 458.9642, 460.8415,
    463.6461, 583.7738, 466.7807, 460.9017, 455.4682, 452.4093, 449.8886, 447.7461,
    446.0029, 444.6366, 443.3456, 442.1168, 441.0445, 439.919, 438.8695, 437.7167,
    436.6012, 435.6349, 434.5094, 433.1968, 431.7222, 430.2362, 428.62, 422.1941,
    282.1162, 256.5904, 240.3565, 233.1687, 227.0982, 224.9568, 219.4109
])

# Step 2: Separate the data into adsorption and desorption branches
# The turning point is where P/P0 is maximum
split_idx = np.argmax(p_p0_raw) + 1
p_ads = p_p0_raw[:split_idx]
v_ads = v_ads_raw[:split_idx]
p_des = p_p0_raw[split_idx:]
v_des = v_ads_raw[split_idx:]

# Step 3: Perform BET analysis to find Vm and SSA
# Filter the adsorption data for the typical BET range (0.05 <= P/P0 <= 0.35)
bet_mask = (p_ads >= 0.05) & (p_ads <= 0.35)
p_bet = p_ads[bet_mask]
v_bet = v_ads[bet_mask]

# Calculate the terms for the BET plot: y = 1 / (V * (P0/P - 1)), x = P/P0
bet_y = 1 / (v_bet * (1/p_bet - 1))
bet_x = p_bet

# Perform linear regression to find slope (m) and intercept (c)
slope, intercept = np.polyfit(bet_x, bet_y, 1)

# Calculate Vm (monolayer capacity) from Vm = 1 / (slope + intercept)
vm = 1 / (slope + intercept)

# Calculate Specific Surface Area (SSA) using the standard conversion factor for N2
ssa = vm * 4.353

# Step 4: Determine the pore diameter (d) from the desorption branch
# Find the largest desorption step (max absolute change in volume)
delta_v = np.abs(np.diff(v_des))
max_des_idx = np.argmax(delta_v)

# Get the P/P0 values for this step and calculate the midpoint pressure
p_step_start = p_des[max_des_idx]
p_step_end = p_des[max_des_idx + 1]
p_midpoint = (p_step_start + p_step_end) / 2

# Calculate the Kelvin radius (r_k) using the Kelvin equation for N2 at 77K
r_k = -0.958 / np.log(p_midpoint)

# Calculate the adsorbed film thickness (t) using the Harkins-Jura equation
# Convert result from Angstroms to nm by dividing by 10
t = 0.1 * (13.99 / (0.034 - np.log10(p_midpoint)))**0.5

# Calculate the pore diameter d = 2 * (r_k + t)
d_calculated = 2 * (r_k + t)

# Round d to the nearest nanometer as requested
d_rounded = round(d_calculated)

# Step 5: Report the results as an ordered set (SSA, Vm, d)
final_answer = (ssa, vm, d_rounded)
print(final_answer)