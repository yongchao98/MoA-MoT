import sys

def analyze_vaccine_efficacy():
    """
    Demonstrates whether 1-IRR correctly estimates per-exposure VE
    for an all-or-nothing vaccine.
    """
    # Step 1: Define a hypothetical scenario
    print("Step 1: Setting up a hypothetical scenario for an 'all-or-nothing' vaccine.")
    population_size = 100000  # Size of each group
    follow_up_years = 1.0     # Follow-up time in years
    # This is the true, underlying efficacy. 80% of vaccinated people are fully immune.
    true_efficacy_aon = 0.80
    # Hazard rate for unprotected individuals (infections per person-year)
    hazard_rate = 0.05

    print(f"- Population per group: {population_size}")
    print(f"- True All-or-Nothing Efficacy (VE_AON): {true_efficacy_aon:.2f} ({true_efficacy_aon:.0%})")
    print(f"- Infection rate for unprotected individuals: {hazard_rate} per person-year\n")

    # Step 2: Calculate expected outcomes and Incidence Rates (IR)

    # --- Unvaccinated Group ---
    person_years_unvaccinated = population_size * follow_up_years
    # For unvaccinated, all are susceptible to the baseline hazard rate.
    cases_unvaccinated = person_years_unvaccinated * hazard_rate
    ir_unvaccinated = cases_unvaccinated / person_years_unvaccinated

    # --- Vaccinated Group ---
    person_years_vaccinated = population_size * follow_up_years
    # In the vaccinated group, only the non-protected fraction is at risk.
    unprotected_fraction = 1 - true_efficacy_aon
    # Cases only occur in the fraction of the vaccinated population that is not protected.
    cases_vaccinated = (population_size * unprotected_fraction) * follow_up_years * hazard_rate
    ir_vaccinated = cases_vaccinated / person_years_vaccinated

    print("Step 2 & 3: Calculate Vaccine Efficacy using the Incidence Rate Ratio (IRR).")
    print("IR = Total Cases / Total Person-Time")
    print(f"- Incidence Rate (Unvaccinated): {int(cases_unvaccinated)} cases / {int(person_years_unvaccinated)} person-years = {ir_unvaccinated:.4f}")
    print(f"- Incidence Rate (Vaccinated):   {int(cases_vaccinated)} cases / {int(person_years_vaccinated)} person-years = {ir_vaccinated:.4f}")

    # Step 3: Calculate VE from 1 - IRR
    irr = ir_vaccinated / ir_unvaccinated
    ve_from_irr = 1 - irr
    
    print("\nThe Vaccine Efficacy calculated from the IRR is:")
    # The format below fulfills the requirement: "output each number in the final equation"
    print(f"VE_IRR = 1 - (IR_Vaccinated / IR_Unvaccinated)")
    sys.stdout.write(f"VE_IRR = 1 - ({ir_vaccinated:.4f} / {ir_unvaccinated:.4f})")
    sys.stdout.flush()
    print(f" = 1 - {irr:.2f} = {ve_from_irr:.2f} ({ve_from_irr:.0%})\n")


    # Step 4: Define and calculate the "true" per-exposure efficacy
    print("Step 4: Calculate the true 'per-exposure' Vaccine Efficacy (VE_p).")
    print("VE_p is the average reduction in infection risk per exposure for a vaccinated person.")
    print("- For the protected fraction (80%), the risk reduction is 100%.")
    print("- For the unprotected fraction (20%), the risk reduction is 0%.")
    
    # Let p be the risk per exposure for an unvaccinated person.
    # The average risk per exposure for a vaccinated person is:
    # (risk_for_protected * fraction_protected) + (risk_for_unprotected * fraction_unprotected)
    # = (0 * true_efficacy_aon) + (p * (1 - true_efficacy_aon)) = p * (1 - true_efficacy_aon)
    # The per-exposure efficacy is 1 - (avg_risk_vaccinated / risk_unvaccinated)
    # VE_p = 1 - (p * (1 - true_efficacy_aon) / p) = 1 - (1 - true_efficacy_aon) = true_efficacy_aon
    
    print("\nThe per-exposure efficacy is mathematically defined as:")
    print("VE_p = 1 - (Average_Risk_Vaccinated / Risk_Unvaccinated)")
    sys.stdout.write(f"VE_p = 1 - (1 - VE_AON)")
    sys.stdout.flush()
    print(f" = 1 - (1 - {true_efficacy_aon:.2f}) = {true_efficacy_aon:.2f} ({true_efficacy_aon:.0%})\n")

    # Step 5: Compare and conclude
    print("Step 5: Conclusion.")
    print(f"The calculated VE from rates (VE_IRR = {ve_from_irr:.2f}) is equal to the true per-exposure efficacy (VE_p = {true_efficacy_aon:.2f}).")
    print("Therefore, for an all-or-nothing vaccine, 1-IRR correctly estimates the per-exposure vaccine efficacy.")


analyze_vaccine_efficacy()
<<<C>>>