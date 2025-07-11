import textwrap

def analyze_all_or_nothing_vaccine():
    """
    Analyzes whether 1 - IRR correctly estimates the efficacy
    of an all-or-nothing vaccine by demonstrating with a hypothetical scenario.
    """
    # Step 1: Define parameters for the hypothetical scenario.
    # True efficacy (VE_s): The proportion of vaccinated people who are fully protected.
    ve_s_true = 0.90  # 90% of vaccinated individuals are rendered fully immune.

    # Force of infection (lambda): The rate at which susceptible individuals get infected.
    # We'll use a rate of 50 cases per 1,000 person-years.
    force_of_infection = 0.050

    # Total person-time observed in each group (e.g., 100,000 people followed for 1 year).
    person_years_unvaccinated = 100000.0
    person_years_vaccinated = 100000.0

    print("--- Scenario: All-or-Nothing Vaccine Efficacy Analysis ---")
    print(f"True Vaccine Efficacy (VE_s), i.e., the proportion protected: {ve_s_true:.2f}")
    print(f"Infection rate for any susceptible person (lambda): {force_of_infection:.3f} cases per person-year")
    print("-" * 60)

    # Step 2: Calculate the incidence rate in the unvaccinated group.
    # In the unvaccinated group, all 100,000 person-years are at risk.
    cases_unvaccinated = person_years_unvaccinated * force_of_infection
    ir_unvaccinated = cases_unvaccinated / person_years_unvaccinated
    print("1. Unvaccinated Group Analysis:")
    print(f"   - Expected Cases = {person_years_unvaccinated:,.0f} person-years * {force_of_infection:.3f} = {cases_unvaccinated:,.0f}")
    print(f"   - Incidence Rate (IR_unvaccinated) = {cases_unvaccinated:,.0f} / {person_years_unvaccinated:,.0f} = {ir_unvaccinated:.3f}")
    print()

    # Step 3: Calculate the incidence rate in the vaccinated group.
    # For an all-or-nothing vaccine, a proportion (1 - VE_s) remains susceptible.
    proportion_susceptible_in_vaccinated = 1 - ve_s_true
    
    # Cases only occur in the susceptible portion of the vaccinated group.
    cases_vaccinated = (person_years_vaccinated * proportion_susceptible_in_vaccinated) * force_of_infection
    
    # The incidence rate is calculated over the *entire* vaccinated group's person-time.
    ir_vaccinated = cases_vaccinated / person_years_vaccinated
    print("2. Vaccinated Group Analysis:")
    print(f"   - Proportion of vaccinated group not protected = 1 - {ve_s_true:.2f} = {proportion_susceptible_in_vaccinated:.2f}")
    print(f"   - Expected Cases = ({person_years_vaccinated:,.0f} * {proportion_susceptible_in_vaccinated:.2f}) * {force_of_infection:.3f} = {cases_vaccinated:,.0f}")
    print(f"   - Incidence Rate (IR_vaccinated) = {cases_vaccinated:,.0f} / {person_years_vaccinated:,.0f} = {ir_vaccinated:.4f}")
    print()

    # Step 4: Calculate the Incidence Rate Ratio (IRR).
    irr = ir_vaccinated / ir_unvaccinated
    print("3. Calculate Incidence Rate Ratio (IRR):")
    print(f"   - IRR = IR_vaccinated / IR_unvaccinated")
    print(f"   - IRR = {ir_vaccinated:.4f} / {ir_unvaccinated:.3f} = {irr:.2f}")
    print()

    # Step 5: Calculate the estimated VE from the IRR.
    ve_estimated = 1 - irr
    print("4. Estimate Vaccine Efficacy (VE) using the 1 - IRR formula:")
    print(f"   - VE_estimated = 1 - IRR")
    print(f"   - VE_estimated = 1 - {irr:.2f} = {ve_estimated:.2f}")
    print()

    # Step 6: Compare the estimated VE to the true VE.
    print("5. Conclusion:")
    print(f"   - The true vaccine efficacy (VE_s) was {ve_s_true:.2f}.")
    print(f"   - The efficacy estimated using the formula 1 - IRR is {ve_estimated:.2f}.")
    
    conclusion_text = ""
    if round(ve_estimated, 4) == round(ve_s_true, 4):
        conclusion_text = "The estimated efficacy is equal to the true efficacy. Therefore, for an all-or-nothing vaccine, 1 - IRR correctly estimates the vaccine efficacy."
    elif ve_estimated > ve_s_true:
        conclusion_text = "The estimated efficacy is greater than the true efficacy. Therefore, for an all-or-nothing vaccine, 1 - IRR overestimates the vaccine efficacy."
    else:
        conclusion_text = "The estimated efficacy is less than the true efficacy. Therefore, for an all-or-nothing vaccine, 1 - IRR underestimates the vaccine efficacy."

    wrapper = textwrap.TextWrapper(width=60, initial_indent='   ', subsequent_indent='   ')
    print(wrapper.fill(conclusion_text))

# Execute the analysis function
analyze_all_or_nothing_vaccine()
<<<C>>>