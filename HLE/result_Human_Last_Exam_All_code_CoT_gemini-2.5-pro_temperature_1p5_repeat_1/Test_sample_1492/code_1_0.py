import math

def solve_shakshuka_dilemma():
    """
    Calculates the required volume of a new shakshuka sauce mixture
    based on changes in ingredient properties.
    """

    # --- Define constants based on the problem statement ---
    
    # Original recipe values
    original_egg_displacement_vol = 45  # in ml (since 1 cm³ = 1 ml)
    
    # New condition values
    num_eggs = 6
    egg_size_increase_percent = 12
    new_consistency_volume_baseline = 180  # in ml

    # --- Perform calculations step-by-step ---

    # Step 1: Calculate the displacement volume for one new, larger egg.
    # The new eggs are 12% larger, so their displacement increases by the same factor.
    egg_size_increase_factor = 1 + (egg_size_increase_percent / 100)
    new_single_egg_displacement = original_egg_displacement_vol * egg_size_increase_factor

    # Step 2: Calculate the total displacement volume for all eggs.
    total_new_displacement_volume = num_eggs * new_single_egg_displacement

    # Step 3: The new volume for sauce consistency is given directly by the 180ml baseline.
    # This represents the "perfect portion" from which the total volume is calculated.

    # Step 4: Sum the new displacement and new consistency volumes for the final answer.
    final_total_volume = total_new_displacement_volume + new_consistency_volume_baseline

    # --- Print the final equation as requested ---
    
    print("To find the total new sauce volume, we sum the volume displaced by the eggs and the volume required for consistency.")
    print("\nCalculating the Final Equation:")
    
    # The final equation is: Total Volume = (Number of Eggs * New Egg Displacement) + New Consistency Volume
    # We will display the calculated values for each part.
    
    final_equation_str = (
        f"Final Volume = ({num_eggs} * {original_egg_displacement_vol} * {egg_size_increase_factor}) + {new_consistency_volume_baseline}"
    )
    
    final_calculation_str = (
        f"Final Volume = {total_new_displacement_volume:.1f} ml + {new_consistency_volume_baseline} ml = {final_total_volume:.1f} ml"
    )

    print(final_equation_str)
    print(final_calculation_str)

solve_shakshuka_dilemma()
<<<482.4>>>