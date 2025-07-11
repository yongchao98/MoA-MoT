def analyze_culture_outcome(media_freshness_factor, incubation_days, num_plates):
    """
    Analyzes the competition between a pathogen and a contaminant on a selective medium.

    Args:
      media_freshness_factor (float): A number from 0 (old/ineffective) to 1 (fresh/perfect).
      incubation_days (int): The number of days the plates are incubated.
      num_plates (int): The number of plates used.
    """
    # Assign base growth scores. The contaminant is naturally a faster grower.
    pathogen_base_growth = 2.0  # e.g., Campylobacter
    contaminant_base_growth = 4.0 # e.g., Bacillus

    # Fresh media suppresses the contaminant and helps the pathogen.
    # We model this by adjusting the growth scores.
    pathogen_final_growth = pathogen_base_growth * media_freshness_factor
    # The contaminant's growth is inversely proportional to media freshness.
    contaminant_final_growth = contaminant_base_growth * (1 - media_freshness_factor)

    # Calculate total colony presence after incubation.
    pathogen_total = pathogen_final_growth * incubation_days * num_plates
    contaminant_total = contaminant_final_growth * incubation_days * num_plates

    print("Analysis based on provided factors:")
    print("Pathogen Score Equation:")
    print(f"    ({pathogen_base_growth} [base] * {media_freshness_factor} [freshness]) * {incubation_days} [days] * {num_plates} [plates] = {pathogen_total:.2f}")
    print("Contaminant Score Equation:")
    print(f"    ({contaminant_base_growth} [base] * (1 - {media_freshness_factor}) [freshness]) * {incubation_days} [days] * {num_plates} [plates] = {contaminant_total:.2f}")

    if pathogen_total > contaminant_total:
        print("\nOutcome: SUCCESS. The pathogen was isolated.")
    else:
        print("\nOutcome: FAILURE. The contaminant outgrew the pathogen.")
    print("-" * 50)


# Scenario 1: The original lab's result, implying old plates with poor selectivity.
print("Original Lab Scenario (assumed old plates):")
# Let's assume old plates have a freshness factor of 0.3 (30% effective).
analyze_culture_outcome(media_freshness_factor=0.3, incubation_days=2, num_plates=1)

# Scenario 2: Applying Solution C - Using fresh plates.
print("Proposed Solution 'C' (using fresh plates):")
# Fresh plates have a high freshness factor, e.g., 0.9 (90% effective).
analyze_culture_outcome(media_freshness_factor=0.9, incubation_days=2, num_plates=1)

<<<C>>>