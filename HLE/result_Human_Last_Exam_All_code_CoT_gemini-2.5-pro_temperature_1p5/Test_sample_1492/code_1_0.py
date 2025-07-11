import math

# Step 1: Define initial state variables from the problem description.
v_old = 750.0  # ml, original total sauce volume
v_disp_per_egg_old = 45.0  # ml, original displacement per egg
num_eggs = 6
v_perfect_portion = 180.0 # ml, the baseline for logarithmic calculation

# Calculate total original egg displacement
v_disp_old = num_eggs * v_disp_per_egg_old

# Step 2: Calculate the change in volume due to larger eggs.
egg_size_increase_factor = 1.12 # 12% larger
v_disp_new = v_disp_old * egg_size_increase_factor
delta_v_disp = v_disp_new - v_disp_old

# Step 3: Calculate the volume reduction due to increased sauce viscosity.
# Viscosity increases by 3/7, so the ratio of new to old is 1 + 3/7 = 10/7.
viscosity_ratio_num = 10
viscosity_ratio_den = 7
viscosity_ratio = viscosity_ratio_num / viscosity_ratio_den

# The reduction is logarithmic, using the 180ml portion as the baseline.
# We use the natural logarithm (ln) as is standard.
delta_v_reduct = v_perfect_portion * math.log(viscosity_ratio)

# Step 4: Calculate the final new sauce volume.
# V_new = V_old + volume_increase_from_eggs - volume_decrease_from_efficiency
v_new = v_old + delta_v_disp - delta_v_reduct

# Step 5: Print the explanation and the final equation with all numbers.
print("To find the new required sauce volume, we adjust the original volume based on two factors:")
print("1. The increase in volume needed to displace the larger eggs.")
print("2. The decrease in volume due to the higher efficiency of the new, more viscous sauce.")
print("\nThe final calculation is:")
print(f"New Volume = Original Volume + Egg Displacement Change - Viscosity Efficiency Savings")
print(f"New Volume = {v_old} ml + ({v_disp_new:.1f} ml - {v_disp_old:.1f} ml) - {v_perfect_portion:.1f} ml * ln({viscosity_ratio_num}/{viscosity_ratio_den})")
print(f"New Volume = {v_old} ml + {delta_v_disp:.1f} ml - {delta_v_reduct:.1f} ml")
print(f"New Volume = {v_new:.1f} ml")

# The final numerical answer requested by the user prompt.
print(f"\n<<<Chef Sarah must prepare {v_new:.1f} ml of the new sauce mixture.>>>")