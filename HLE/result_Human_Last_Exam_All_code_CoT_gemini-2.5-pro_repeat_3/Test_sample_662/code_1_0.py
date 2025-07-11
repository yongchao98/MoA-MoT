import sys

def calculate_vaccine_coverage():
    """
    Calculates the theoretical threshold of vaccine coverage needed to prevent an epidemic.
    """
    # Known values from the user's problem description
    R0 = 3.0  # Basic Reproduction Number
    infection_rate_in_vaccinated_percent = 6.0

    # Convert percentage to a decimal for calculation
    infection_rate_in_vaccinated = infection_rate_in_vaccinated_percent / 100.0

    # Step 1: Calculate Vaccine Effectiveness (VE)
    # If 6% of vaccinated people get infected, the vaccine is effective in the other 94%.
    vaccine_effectiveness = 1 - infection_rate_in_vaccinated

    # Step 2: Calculate the Herd Immunity Threshold (HIT)
    # This is the proportion of the population that needs to be immune.
    # Formula: HIT = 1 - (1 / R0)
    herd_immunity_threshold = 1 - (1 / R0)

    # Step 3: Calculate the Critical Vaccination Coverage (p_c)
    # This adjusts the HIT for the vaccine's effectiveness.
    # Formula: p_c = HIT / VE
    critical_coverage = herd_immunity_threshold / vaccine_effectiveness

    # --- Output the results ---
    print("Calculating the Theoretical Vaccine Coverage Threshold")
    print("-" * 50)
    print(f"Given values:")
    print(f"  - Basic Reproduction Number (R₀): {R0}")
    print(f"  - Infection Rate in Vaccinated People: {infection_rate_in_vaccinated_percent}%")
    print("\nStep 1: Determine Vaccine Effectiveness (VE)")
    print(f"   VE = 100% - (Infection Rate in Vaccinated)")
    print(f"   VE = 1 - {infection_rate_in_vaccinated} = {vaccine_effectiveness:.2f} or {vaccine_effectiveness:.0%}")
    
    print("\nStep 2: Determine Herd Immunity Threshold (HIT)")
    print(f"   HIT = 1 - (1 / R₀)")
    print(f"   HIT = 1 - (1 / {R0}) = {herd_immunity_threshold:.4f} or {herd_immunity_threshold:.1%}")

    print("\nStep 3: Calculate Critical Vaccination Coverage (p_c)")
    print("   This is the coverage needed to achieve herd immunity with a less-than-perfect vaccine.")
    print(f"   p_c = HIT / VE")
    
    # Final equation with numbers, as requested
    print("\nFinal Equation:")
    print(f"   Required Coverage = (1 - (1 / {R0})) / {vaccine_effectiveness}")
    print(f"   Required Coverage = {herd_immunity_threshold:.4f} / {vaccine_effectiveness:.2f}")
    print(f"   Required Coverage = {critical_coverage:.4f}")

    print("\n" + "=" * 50)
    print(f"Conclusion: To prevent the spread of the virus, the theoretical vaccine coverage threshold is {critical_coverage:.1%}.")
    print("=" * 50)

    # Final answer in the required format for automated checking.
    final_answer = critical_coverage * 100
    sys.stdout.write(f"\n<<<{final_answer:.1f}>>>\n")

if __name__ == '__main__':
    calculate_vaccine_coverage()