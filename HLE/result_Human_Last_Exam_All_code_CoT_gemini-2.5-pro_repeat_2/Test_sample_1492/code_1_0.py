import math

def solve_shakshuka_dynamics():
    """
    Calculates the new sauce volume for Chef Sarah's shakshuka based on changes
    in ingredient properties.
    """
    # Step 1: Define initial recipe components
    original_total_volume = 750.0  # ml
    num_eggs = 6
    original_displacement_per_egg = 45.0  # ml
    perfect_portion_baseline = 180.0  # ml

    # Step 2: Deconstruct the original sauce volume
    original_total_displacement = num_eggs * original_displacement_per_egg
    original_free_sauce = original_total_volume - original_total_displacement
    
    # The free sauce is composed of the baseline and an adjustable portion
    original_adjustable_sauce = original_free_sauce - perfect_portion_baseline

    # Step 3: Calculate the new volume required for egg displacement
    egg_size_increase_factor = 1.12  # 12% larger
    new_total_displacement = original_total_displacement * egg_size_increase_factor

    # Step 4: Calculate the new volume for the adjustable portion of the free sauce
    # Viscosity increases by 3/7, so the ratio of new to old is (1 + 3/7) = 10/7.
    # The volume requirement changes logarithmically with this ratio.
    viscosity_ratio = 10.0 / 7.0
    logarithmic_reduction_factor = math.log(viscosity_ratio)
    
    # The reduction applies to the original adjustable sauce volume
    volume_reduction = original_adjustable_sauce * logarithmic_reduction_factor
    new_adjustable_sauce = original_adjustable_sauce - volume_reduction
    
    # Step 5: Calculate the final total sauce volume
    # Final Volume = Baseline + New Adjustable Sauce + New Displacement
    final_total_volume = perfect_portion_baseline + new_adjustable_sauce + new_total_displacement

    # Output the final equation and result
    print("Calculating the new required sauce volume:")
    print(f"New Volume = (Perfect Portion) + (Original Adjustable Sauce - Reduction) + (New Egg Displacement)")
    print(f"New Volume = {perfect_portion_baseline:.1f} ml + ({original_adjustable_sauce:.1f} ml - {volume_reduction:.1f} ml) + ({new_total_displacement:.1f} ml)")
    print(f"New Volume = {perfect_portion_baseline:.1f} ml + {new_adjustable_sauce:.1f} ml + {new_total_displacement:.1f} ml")
    print(f"Final Required Sauce Volume = {final_total_volume:.1f} ml")

solve_shakshuka_dynamics()