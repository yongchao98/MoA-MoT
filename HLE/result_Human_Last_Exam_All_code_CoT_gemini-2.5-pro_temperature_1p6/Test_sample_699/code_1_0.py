import math

def estimate_fatigue_life():
    """
    Estimates fatigue life based on the Basquin relationship and Palmgren-Miner rule.
    """
    # --- Step 1: Define constants and the key assumption ---
    b = -0.09  # Basquin exponent

    # This is a key assumption: The endurance limit is half the ultimate tensile strength.
    # This is required to get a numerical answer as the ratio was not provided.
    uts_to_se_ratio = 0.5

    # Loading profile: cycle percentages and corresponding stress multipliers
    loading_profile = {
        "level_1": {"cycle_percent": 0.70, "stress_multiplier": 1.0},
        "level_2": {"cycle_percent": 0.20, "stress_multiplier": 1.1},
        "level_3": {"cycle_percent": 0.10, "stress_multiplier": 1.2},
    }

    # --- Step 2: Calculate N_f (cycles to failure) for each stress level ---
    # The general S-N curve equation derived from the problem statement is:
    # N_f = 0.5 * (sigma_a / UTS)^(1/b)
    # Since sigma_a = stress_multiplier * sigma_e and sigma_e = uts_to_se_ratio * UTS,
    # the ratio (sigma_a / UTS) becomes (stress_multiplier * uts_to_se_ratio).

    cycles_to_failure = {}
    for level, data in loading_profile.items():
        stress_ratio_vs_uts = data["stress_multiplier"] * uts_to_se_ratio
        N_f = 0.5 * (stress_ratio_vs_uts)**(1 / b)
        cycles_to_failure[level] = N_f

    N_f1 = cycles_to_failure["level_1"]
    N_f2 = cycles_to_failure["level_2"]
    N_f3 = cycles_to_failure["level_3"]
    
    c1 = loading_profile["level_1"]["cycle_percent"]
    c2 = loading_profile["level_2"]["cycle_percent"]
    c3 = loading_profile["level_3"]["cycle_percent"]

    # --- Step 3: Apply Palmgren-Miner Rule ---
    # The total damage from one repeating block of cycles is the sum of partial damages.
    # N_life = 1 / Sum(cycle_percent_i / N_fi)
    
    damage_sum = (c1 / N_f1) + (c2 / N_f2) + (c3 / N_f3)
    
    if damage_sum <= 0:
        N_life = float('inf')
    else:
        N_life = 1 / damage_sum

    # --- Step 4: Print the full calculation and the final answer ---
    print("--- Fatigue Life Estimation ---")
    print(f"\nGiven Basquin exponent (b): {b}")
    print(f"Assumed Endurance Limit (σ_e) = {uts_to_se_ratio} * UTS")

    print("\nCycles to Failure at each stress level:")
    print(f"N_f1 (at 1.0σ_e): {N_f1:.3e}")
    print(f"N_f2 (at 1.1σ_e): {N_f2:.3e}")
    print(f"N_f3 (at 1.2σ_e): {N_f3:.3e}")

    print("\nApplying the Palmgren-Miner rule:")
    print("N_life = 1 / [ (% of cycles at level 1 / N_f1) + (% of cycles at level 2 / N_f2) + (% of cycles at level 3 / N_f3) ]")
    print(f"N_life = 1 / [ ({c1} / {N_f1:.3e}) + ({c2} / {N_f2:.3e}) + ({c3} / {N_f3:.3e}) ]")
    print(f"N_life = 1 / [ {c1/N_f1:.3e} + {c2/N_f2:.3e} + {c3/N_f3:.3e} ]")
    print(f"N_life = 1 / [ {damage_sum:.3e} ]")
    
    print("\n------------------------------------")
    print("Final Estimated Fatigue Life:")
    print(f"{N_life:.3e} cycles")
    print("------------------------------------")
    
    # Return the formatted string for the final answer block
    return f"{N_life:.3e}"

# Execute the function and capture the final answer
final_answer = estimate_fatigue_life()
print(f"<<<{final_answer}>>>")
