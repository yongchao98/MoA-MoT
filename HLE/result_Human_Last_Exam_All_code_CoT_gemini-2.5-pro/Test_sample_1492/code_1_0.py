import math

# Step 1: Define the initial values from the problem description.
original_sauce_volume_ml = 750
perfect_portion_ml = 180
viscosity_increase_fraction_num = 3
viscosity_increase_fraction_den = 7
num_eggs = 6
original_egg_displacement_ml = 45
egg_size_increase = 0.12 # 12%

# Step 2: Calculate the portion of the original sauce volume affected by the changes.
# This is the original total volume minus the "perfect" portion.
affected_volume_ml = original_sauce_volume_ml - perfect_portion_ml

# Step 3: Calculate the adjusted volume for the affected portion.
# The viscosity increased by 3/7, leading to a proportional efficiency gain.
# This means the required volume for this portion is reduced by 3/7.
efficiency_multiplier = 1 - (viscosity_increase_fraction_num / viscosity_increase_fraction_den)
adjusted_affected_volume_ml = affected_volume_ml * efficiency_multiplier

# Step 4: Calculate the additional volume required due to larger eggs.
# The new eggs are 12% larger, so they displace more sauce.
original_total_displacement = num_eggs * original_egg_displacement_ml
new_total_displacement = original_total_displacement * (1 + egg_size_increase)
extra_volume_for_eggs_ml = new_total_displacement - original_total_displacement

# Step 5: Sum the components to find the final required sauce volume.
# Final Volume = (Unaffected Portion) + (Adjusted Affected Portion) + (Extra Egg Displacement)
final_volume_ml = perfect_portion_ml + adjusted_affected_volume_ml + extra_volume_for_eggs_ml

# Step 6: Print the final equation with all the original numbers and the result.
# The request requires showing each number in the final equation.
print("To calculate the required new sauce volume, we apply adjustments to the original 750ml recipe:")
print("Final Volume = (Perfect Portion) + (Original Affected Portion * Efficiency Factor) + (Extra Egg Displacement)")
print(f"Final Equation: {perfect_portion_ml} + ({original_sauce_volume_ml} - {perfect_portion_ml}) * (1 - {viscosity_increase_fraction_num}/{viscosity_increase_fraction_den}) + (({num_eggs} * {original_egg_displacement_ml} * (1 + {egg_size_increase})) - ({num_eggs} * {original_egg_displacement_ml})) = {final_volume_ml:.1f} ml")

print(f"\nChef Sarah must prepare approximately {round(final_volume_ml)} ml of the new sauce.")

# Final answer in the required format.
# The calculation is 180 + (570 * 4/7) + 32.4 = 180 + 325.71... + 32.4 = 538.11...
final_answer = round(final_volume_ml)
# The final answer is wrapped in <<<>>>
# <<<538>>>