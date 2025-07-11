import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    with common cause failures to meet an overall ASIL C target.
    """
    # --- Given Parameters ---
    lambda_ods1_fit = 500  # Failure rate of ODS1 in FIT
    lambda_ods2_fit = 400  # Failure rate of ODS2 in FIT
    lambda_ods3_fit = 700  # Failure rate of ODS3 in FIT
    beta = 0.10            # Common Cause Failure factor
    t_mission_h = 10000    # Mission time in hours
    lambda_target_fit = 100 # System failure rate target for ASIL C in FIT

    # Conversion factor from FIT to failures per hour
    fit_to_per_hour = 1e-9

    print("Step 1: Calculate the Common Cause Failure (CCF) contribution (λ_ccf).")
    lambda_ods_sum_fit = lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit
    lambda_ods_avg_fit = lambda_ods_sum_fit / 3
    lambda_ccf_fit = beta * lambda_ods_avg_fit
    print(f"λ_ccf = β * (λ_ODS1 + λ_ODS2 + λ_ODS3) / 3")
    print(f"λ_ccf = {beta} * ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) / 3 = {lambda_ccf_fit:.3f} FIT\n")

    print("Step 2: Calculate the contribution from independent failures (λ_2oo3_ind).")
    # Independent failure rates in FIT
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit
    print(f"λ_ODS1_ind = (1 - {beta}) * {lambda_ods1_fit} = {lambda_ods1_ind_fit} FIT")
    print(f"λ_ODS2_ind = (1 - {beta}) * {lambda_ods2_fit} = {lambda_ods2_ind_fit} FIT")
    print(f"λ_ODS3_ind = (1 - {beta}) * {lambda_ods3_fit} = {lambda_ods3_ind_fit} FIT\n")

    # Convert independent rates from FIT to per-hour for the formula
    lambda_ods1_ind_hr = lambda_ods1_ind_fit * fit_to_per_hour
    lambda_ods2_ind_hr = lambda_ods2_ind_fit * fit_to_per_hour
    lambda_ods3_ind_hr = lambda_ods3_ind_fit * fit_to_per_hour

    # Sum of the products of paired independent failure rates
    pair_product_sum_hr2 = (lambda_ods1_ind_hr * lambda_ods2_ind_hr +
                            lambda_ods1_ind_hr * lambda_ods3_ind_hr +
                            lambda_ods2_ind_hr * lambda_ods3_ind_hr)

    # Average failure rate of the 2oo3 independent system over the mission time
    lambda_2oo3_ind_avg_hr = pair_product_sum_hr2 * t_mission_h

    # Convert the result back to FIT
    lambda_2oo3_ind_avg_fit = lambda_2oo3_ind_avg_hr / fit_to_per_hour
    print(f"Calculating the average failure rate from independent pairs over {t_mission_h} hours:")
    print(f"λ_2oo3_ind ≈ (λ_1_ind*λ_2_ind + λ_1_ind*λ_3_ind + λ_2_ind*λ_3_ind) * T_mission")
    print(f"λ_2oo3_ind = {lambda_2oo3_ind_avg_fit:.3f} FIT\n")

    print("Step 3: Determine the voter failure rate budget (λ_voter).")
    # Calculate the remaining budget for the voter
    lambda_voter_budget_fit = lambda_target_fit - lambda_ccf_fit - lambda_2oo3_ind_avg_fit
    
    print("The total system failure rate must be less than the target.")
    print("λ_system = λ_2oo3_ind + λ_ccf + λ_voter < λ_target")
    print(f"Substituting the calculated values:")
    print(f"{lambda_2oo3_ind_avg_fit:.3f} FIT + {lambda_ccf_fit:.3f} FIT + λ_voter < {lambda_target_fit} FIT")
    print(f"Solving for λ_voter:")
    print(f"λ_voter < {lambda_target_fit} - {lambda_ccf_fit:.3f} - {lambda_2oo3_ind_avg_fit:.3f}")
    print(f"λ_voter < {lambda_voter_budget_fit:.3f} FIT\n")

    # Final Answer
    print("The required failure rate for the voter is:")
    print(f"λvoter < {lambda_voter_budget_fit:.3f} FIT")
    
    return lambda_voter_budget_fit

# Execute the calculation
final_answer = calculate_voter_failure_rate()
print(f"<<<{final_answer:.3f}>>>")
