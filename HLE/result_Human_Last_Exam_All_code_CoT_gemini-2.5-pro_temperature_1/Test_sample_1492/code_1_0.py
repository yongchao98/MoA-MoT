import math

# Step 1: Define initial variables from the problem description
v_total_orig = 750  # Original total sauce volume in ml
v_baseline = 180    # Baseline volume that serves as the mathematical reference
egg_increase = 0.12 # The percentage increase in egg size
viscosity_increase_fraction_numerator = 3
viscosity_increase_fraction_denominator = 7

print("Solving Chef Sarah's Shakshuka Dilemma\n")
print(f"Original total sauce volume: {v_total_orig} ml")
print(f"Baseline sauce volume: {v_baseline} ml")

# Step 2: Calculate the initial adjustable volume
v_adj_orig = v_total_orig - v_baseline
print(f"Initial adjustable volume = {v_total_orig} - {v_baseline} = {v_adj_orig} ml")

# Step 3: Adjust the volume for the larger eggs
egg_multiplier = 1 + egg_increase
v_adj_for_eggs = v_adj_orig * egg_multiplier
print(f"\nVolume adjusted for 12% larger eggs = {v_adj_orig} * {egg_multiplier} = {v_adj_for_eggs:.2f} ml")

# Step 4: Adjust the volume for the change in sauce viscosity
# The viscosity increases by 3/7, so the new viscosity is 1 + 3/7 = 10/7 times the old viscosity.
viscosity_ratio = 1 + (viscosity_increase_fraction_numerator / viscosity_increase_fraction_denominator)

# The problem states the volume requirement changes logarithmically with viscosity.
# Increased viscosity means higher efficiency, so less volume is needed.
# We model this as a reduction factor based on the natural log of the viscosity ratio.
log_reduction = math.log(viscosity_ratio)
efficiency_multiplier = 1 - log_reduction

print(f"\nViscosity increases by {viscosity_increase_fraction_numerator}/{viscosity_increase_fraction_denominator}, for a ratio of {viscosity_ratio:.2f}.")
print(f"Logarithmic efficiency multiplier = 1 - ln({viscosity_ratio:.2f}) = {efficiency_multiplier:.4f}")

v_adj_final = v_adj_for_eggs * efficiency_multiplier
print(f"Final adjustable volume = {v_adj_for_eggs:.2f} * {efficiency_multiplier:.4f} = {v_adj_final:.2f} ml")

# Step 5: Calculate the final total volume
v_total_new = v_baseline + v_adj_final
print(f"\nFinal total sauce volume = Baseline Volume + Final Adjustable Volume")
print(f"Final total sauce volume = {v_baseline} + {v_adj_final:.2f} = {v_total_new:.2f} ml")

print("\n--- Final Calculation Summary ---")
print("The final volume is calculated using the following equation:")
final_equation = f"Final Volume = {v_baseline} + ({v_total_orig} - {v_baseline}) * (1 + {egg_increase}) * (1 - ln(1 + {viscosity_increase_fraction_numerator}/{viscosity_increase_fraction_denominator}))"
print(final_equation)
print(f"Result: {v_total_new:.2f} ml")

# The final answer rounded to two decimal places
final_answer = round(v_total_new, 2)
# print(f"\n<<<{final_answer}>>>")