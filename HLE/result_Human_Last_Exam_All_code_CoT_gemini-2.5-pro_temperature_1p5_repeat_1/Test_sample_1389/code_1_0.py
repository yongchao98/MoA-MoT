def analyze_all_or_nothing_vaccine():
    """
    Analyzes if 1 - IRR correctly estimates per-exposure efficacy
    for an all-or-nothing vaccine.
    """
    # Step 1: Define the true per-exposure vaccine efficacy (VE_p).
    # For an all-or-nothing vaccine, this is the proportion of people who are
    # fully protected. Let's assume 85% are fully protected.
    ve_p_true_proportion = 0.85
    print(f"Let's assume the true per-exposure Vaccine Efficacy (VE_p) is {ve_p_true_proportion:.2f}.")
    print("This means the vaccine gives perfect protection to 85% of recipients, and no protection to 15%.")
    print("-" * 50)

    # Step 2: Define a baseline incidence rate in the unvaccinated (IR_u).
    # This can be any arbitrary rate, e.g., 10 cases per 1000 person-years.
    # The actual value doesn't matter, as it will cancel out.
    ir_unvaccinated = 0.010
    print(f"Let's assume the incidence rate in the unvaccinated group (IR_u) is {ir_unvaccinated:.3f}.")
    print("-" * 50)

    # Step 3: Calculate the incidence rate in the vaccinated group (IR_v).
    # Since 85% of the vaccinated are perfectly protected, their incidence rate is 0.
    # The remaining 15% (1 - 0.85) are unprotected and have the same incidence
    # rate as the unvaccinated group.
    proportion_unprotected = 1 - ve_p_true_proportion
    ir_vaccinated = ir_unvaccinated * proportion_unprotected
    print("The incidence rate in the vaccinated group (IR_v) is the rate among the unprotected proportion.")
    print(f"IR_v = IR_u * (1 - VE_p)")
    print(f"IR_v = {ir_unvaccinated:.3f} * (1 - {ve_p_true_proportion:.2f})")
    print(f"IR_v = {ir_unvaccinated:.3f} * {proportion_unprotected:.2f}")
    print(f"IR_v = {ir_vaccinated:.5f}")
    print("-" * 50)

    # Step 4: Calculate the Incidence Rate Ratio (IRR).
    # IRR = IR_v / IR_u
    irr = ir_vaccinated / ir_unvaccinated
    print("Next, we calculate the Incidence Rate Ratio (IRR).")
    print("IRR = IR_v / IR_u")
    print(f"IRR = {ir_vaccinated:.5f} / {ir_unvaccinated:.3f}")
    print(f"IRR = {irr:.2f}")
    print("-" * 50)

    # Step 5: Calculate the estimated VE from the IRR.
    # VE_estimated = 1 - IRR
    ve_estimated = 1 - irr
    print("Now, we calculate the Vaccine Efficacy using the standard formula.")
    print("VE_estimated = 1 - IRR")
    print(f"VE_estimated = 1 - {irr:.2f}")
    print(f"VE_estimated = {ve_estimated:.2f}")
    print("-" * 50)

    # Step 6: Compare the estimated VE with the true VE_p.
    print("Finally, we compare the estimated VE with the true per-exposure VE (VE_p).")
    print(f"True VE_p = {ve_p_true_proportion:.2f}")
    print(f"Estimated VE from 1-IRR = {ve_estimated:.2f}")

    if round(ve_p_true_proportion, 5) == round(ve_estimated, 5):
        print("\nConclusion: The value from '1 - IRR' exactly matches the true per-exposure vaccine efficacy.")
        final_answer = "correctly estimate"
    elif ve_estimated > ve_p_true_proportion:
        print("\nConclusion: '1 - IRR' overestimates the true per-exposure vaccine efficacy.")
        final_answer = "overestimate"
    else:
        print("\nConclusion: '1 - IRR' underestimates the true per-exposure vaccine efficacy.")
        final_answer = "underestimate"

    print(f"\nTherefore, for an all-or-nothing vaccine, 1 - incidence rate ratio will {final_answer} the per-exposure vaccine efficacy.")

analyze_all_or_nothing_vaccine()