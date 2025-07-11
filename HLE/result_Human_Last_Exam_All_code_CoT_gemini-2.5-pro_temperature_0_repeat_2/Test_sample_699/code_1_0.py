import math

# 1. Define constants and assumptions
# Basquin exponent
b = -0.09
# Life fractions for each stress level
fraction_1 = 0.70 # at sigma_e
fraction_2 = 0.20 # at 1.1 * sigma_e
fraction_3 = 0.10 # at 1.2 * sigma_e

# Assumption: sigma_e = 0.5 * sigma_uts.
# This allows us to define the damaging stress levels relative to sigma_uts.
stress_ratio_2 = 1.1 * 0.5  # Corresponds to sigma_a = 1.1 * sigma_e
stress_ratio_3 = 1.2 * 0.5  # Corresponds to sigma_a = 1.2 * sigma_e

# 2. Calculate N_f for the damaging stress levels.
# The general formula for N_f, derived from the problem statement, is:
# N_f = ( (sigma_a / sigma_uts) * (0.5)^b )^(1/b)
# We only need to calculate N_f for the stress levels that cause damage (i.e., above sigma_e).
N_f2 = (stress_ratio_2 * (0.5)**b)**(1/b)
N_f3 = (stress_ratio_3 * (0.5)**b)**(1/b)

# 3. Apply the Palmgren-Miner rule to find the total life (N_total).
# The damage from cycles at the endurance limit (sigma_e) is zero because N_f1 is infinite.
# The damage equation is: (fraction_2 / N_f2) + (fraction_3 / N_f3) = 1 / N_total
# Therefore, N_total = 1 / ( (fraction_2 / N_f2) + (fraction_3 / N_f3) )
total_life = 1 / (fraction_2 / N_f2 + fraction_3 / N_f3)

# 4. Output the results, including the final equation with numbers.
print("Based on the Palmgren-Miner rule and the given assumptions:")
print(f"Cycles to failure at 1.1*σ_e (N_f2): {N_f2:.3f}")
print(f"Cycles to failure at 1.2*σ_e (N_f3): {N_f3:.3f}\n")

print("The final equation for the total fatigue life (N_total) is:")
print(f"N_total = 1 / ( ({fraction_2} / {N_f2:.3f}) + ({fraction_3} / {N_f3:.3f}) )")
print(f"N_total = {total_life:.3f} cycles\n")

# Format the final answer into "Aeb" scientific notation with 3 decimal places.
formatted_life = "{:.3e}".format(total_life)
parts = formatted_life.split('e')
mantissa = parts[0]
exponent = int(parts[1])
final_answer_str = f"{mantissa}e{exponent:+d}"

print(f"The estimated fatigue life in the specified format is: {final_answer_str}")

# Final answer block
print(f"<<<{final_answer_str}>>>")