import sys

def analyze_vaccine_efficacy():
    """
    Analyzes if 1 - IRR correctly estimates the efficacy of an all-or-nothing vaccine.
    """
    # Step 1: Define the model parameters with a numerical example.
    # Let's assume the true vaccine efficacy (VE_true) is 85%, meaning it provides
    # complete protection to 85% of vaccinated individuals.
    ve_true = 0.85

    # Let lambda be the force of infection (incidence rate in a fully susceptible population).
    # We can use an arbitrary value, as it will cancel out. Let's use 0.05 for clarity.
    lambda_rate = 0.05

    print("--- Model Definition ---")
    print(f"True Vaccine Efficacy (VE_true), the proportion fully protected: {ve_true}")
    print(f"Assumed incidence rate among susceptibles (lambda): {lambda_rate}")
    print("-" * 25)

    # Step 2: Calculate the theoretical incidence rates for each group.
    # The incidence rate in the unvaccinated group (IR_u) is lambda.
    ir_u = lambda_rate

    # The incidence rate in the vaccinated group (IR_v) is lambda multiplied by the
    # fraction of vaccinated individuals who are NOT protected (1 - VE_true).
    ir_v = lambda_rate * (1 - ve_true)

    print("--- Calculating Incidence Rates ---")
    print(f"Incidence Rate in Unvaccinated (IR_u) = lambda = {ir_u}")
    print(f"Incidence Rate in Vaccinated (IR_v) = lambda * (1 - VE_true) = {lambda_rate} * (1 - {ve_true}) = {ir_v:.4f}")
    print("-" * 25)

    # Step 3: Calculate the Incidence Rate Ratio (IRR) and the estimated VE.
    # The IRR is the ratio of the two incidence rates.
    irr = ir_v / ir_u

    # The estimated VE is calculated as 1 - IRR.
    ve_estimated = 1 - irr

    print("--- Estimating Vaccine Efficacy ---")
    print(f"Incidence Rate Ratio (IRR) = IR_v / IR_u = {ir_v:.4f} / {ir_u} = {irr:.2f}")
    
    # Final equation with numbers printed out
    print("\nThe final efficacy equation is: VE_estimated = 1 - IRR")
    print(f"So, the estimated efficacy is: {ve_estimated:.2f} = 1 - {irr:.2f}")
    print("-" * 25)

    # Step 4: Compare the estimated VE with the true VE.
    print("--- Conclusion ---")
    print(f"The true vaccine efficacy was set to: {ve_true}")
    print(f"The estimated vaccine efficacy (1 - IRR) is: {ve_estimated:.2f}")

    if abs(ve_estimated - ve_true) < 1e-9:
        print("\nResult: 1 - IRR correctly estimates the per-exposure vaccine efficacy.")
        # We use exit codes to return the final answer choice for automated checking.
        # This will not affect manual execution.
        sys.exit(0) # Corresponds to 'C'
    elif ve_estimated > ve_true:
        print("\nResult: 1 - IRR overestimates the per-exposure vaccine efficacy.")
        sys.exit(1) # Corresponds to 'A'
    else:
        print("\nResult: 1 - IRR underestimates the per-exposure vaccine efficacy.")
        sys.exit(2) # Corresponds to 'B'

# Execute the analysis
try:
    analyze_vaccine_efficacy()
except SystemExit as e:
    # This block is to handle the sys.exit call and map it to the final answer.
    # The print statements inside the function have already provided the detailed explanation.
    pass
