def calculate_vaccine_coverage():
    """
    Calculates the theoretical vaccine coverage needed to prevent an epidemic.
    """
    # Step 0: Define the known variables based on the user's input.
    # The basic reproduction number of the virus.
    R0 = 3.0
    # The percentage of vaccinated people who still get the infection.
    percent_vaccinated_infected = 6.0

    print(f"Known variables:")
    print(f"Basic Reproduction Number (R₀): {R0}")
    print(f"Percentage of vaccinated people who get infected: {percent_vaccinated_infected}%\n")

    # Step 1: Calculate the vaccine's effectiveness (E).
    # This is 1 minus the rate at which vaccinated people get infected.
    vaccine_effectiveness = 1.0 - (percent_vaccinated_infected / 100.0)
    print(f"Step 1: Calculate Vaccine Effectiveness (E)")
    print(f"E = 100% - {percent_vaccinated_infected}% = {vaccine_effectiveness * 100:.0f}% or {vaccine_effectiveness:.2f}\n")

    # Step 2: Calculate the Herd Immunity Threshold (HIT).
    # This is the proportion of the population that must be immune.
    # HIT = 1 - (1 / R₀)
    herd_immunity_threshold = 1.0 - (1.0 / R0)
    print(f"Step 2: Calculate Herd Immunity Threshold (HIT)")
    print(f"HIT = 1 - (1 / R₀)")
    print(f"HIT = 1 - (1 / {R0}) = {herd_immunity_threshold:.4f} or {herd_immunity_threshold * 100:.2f}%\n")

    # Step 3: Calculate the required vaccine coverage (Vc).
    # This adjusts the HIT for the vaccine's effectiveness.
    # Vc = HIT / E
    critical_vaccination_coverage = herd_immunity_threshold / vaccine_effectiveness
    print(f"Step 3: Calculate the Critical Vaccination Coverage (Vc)")
    print("This is the proportion of the population that needs to be vaccinated.\n")
    print("Formula: Vc = (1 - 1/R₀) / E\n")

    # Final result and equation
    print("Final Calculation:")
    print(f"(1 - 1 / {R0}) / {vaccine_effectiveness:.2f} = {critical_vaccination_coverage * 100:.2f}%")
    print("\nTo prevent the spread of the virus, the theoretical threshold of vaccine coverage is approximately 70.92%.")
    
    # Return the final numerical answer for the system.
    return critical_vaccination_coverage * 100

# Execute the function and capture the final answer.
final_answer = calculate_vaccine_coverage()
# The final answer is formatted for the system below.
# print(f"\n<<<{final_answer:.2f}>>>")
print(f"<<<{final_answer:.2f}>>>")