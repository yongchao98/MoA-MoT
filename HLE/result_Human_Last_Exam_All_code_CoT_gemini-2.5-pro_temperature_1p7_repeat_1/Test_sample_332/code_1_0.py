def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3
    redundant system with common cause failures to meet a safety target.
    """
    # Define initial parameters from the problem description
    lambda_ods1_fit = 500
    lambda_ods2_fit = 400
    lambda_ods3_fit = 700
    beta = 0.10
    lambda_system_target_fit = 100
    t_mission_h = 10000

    # Conversion factor from failures per hour to FIT (1 FIT = 1e-9 failures/hour)
    fit_conversion_factor = 1e9

    print("This program calculates the required failure rate for a voter in a 2-out-of-3 system.")
    print("The total system failure rate equation is:")
    print("λ_system = λ_2oo3_independent + λ_common_cause + λ_voter\n")
    print(f"The target for the system (ASIL C) is: λ_system < {lambda_system_target_fit} FIT\n")

    # Step 1: Calculate the Common Cause Failure (CCF) rate (λ_ccf)
    lambda_ods_total_fit = lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit
    lambda_ccf_fit = beta * lambda_ods_total_fit

    print("Step 1: Calculate the Common Cause Failure (CCF) rate (λ_ccf).")
    print(f"λ_ccf = β * (λ_ODS1 + λ_ODS2 + λ_ODS3)")
    print(f"λ_ccf = {beta:.2f} * ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit})")
    print(f"λ_ccf = {beta:.2f} * {lambda_ods_total_fit} = {lambda_ccf_fit:.3f} FIT\n")

    # Step 2: Calculate the failure rate of the 2-out-of-3 system from independent failures (λ_2oo3_independent)
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit
    
    # Convert independent FIT rates to failures per hour for reliability calculation
    lambda_ods1_ind_h = lambda_ods1_ind_fit / fit_conversion_factor
    lambda_ods2_ind_h = lambda_ods2_ind_fit / fit_conversion_factor
    lambda_ods3_ind_h = lambda_ods3_ind_fit / fit_conversion_factor

    # Approximate the unreliability F(T) over the mission time
    # F(T) ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * T^2
    term12 = lambda_ods1_ind_h * lambda_ods2_ind_h
    term13 = lambda_ods1_ind_h * lambda_ods3_ind_h
    term23 = lambda_ods2_ind_h * lambda_ods3_ind_h
    f_2oo3_t = (term12 + term13 + term23) * (t_mission_h ** 2)

    # Average failure rate λ_avg = F(T) / T
    lambda_2oo3_ind_avg_h = f_2oo3_t / t_mission_h
    lambda_2oo3_ind_fit = lambda_2oo3_ind_avg_h * fit_conversion_factor

    print("Step 2: Calculate the failure rate of the 2-out-of-3 system from independent failures (λ_2oo3_independent).")
    print("First, calculate the independent failure rates:")
    print(f"λ_1_ind = (1 - {beta:.2f}) * {lambda_ods1_fit} = {lambda_ods1_ind_fit:.3f} FIT")
    print(f"λ_2_ind = (1 - {beta:.2f}) * {lambda_ods2_fit} = {lambda_ods2_ind_fit:.3f} FIT")
    print(f"λ_3_ind = (1 - {beta:.2f}) * {lambda_ods3_fit} = {lambda_ods3_ind_fit:.3f} FIT")
    print(f"\nNext, the average failure rate for the 2oo3 configuration over a {t_mission_h}h mission is calculated.")
    print(f"λ_2oo3_independent ≈ {lambda_2oo3_ind_fit:.3f} FIT\n")

    # Step 3: Calculate the maximum allowed failure rate for the voter
    max_lambda_voter_fit = lambda_system_target_fit - lambda_ccf_fit - lambda_2oo3_ind_fit

    print("Step 3: Calculate the required failure rate for the voter (λ_voter).")
    print("This is found by rearranging the system failure rate equation:")
    print("λ_voter < λ_system_target - λ_ccf - λ_2oo3_independent")
    print(f"λ_voter < {lambda_system_target_fit} - {lambda_ccf_fit:.3f} - {lambda_2oo3_ind_fit:.3f}")
    print(f"\nFinal Result: λvoter < {max_lambda_voter_fit:.3f} FIT\n")

    print("Note: A negative failure rate is physically impossible. This result indicates")
    print("that the system design cannot meet the ASIL C target of 100 FIT, primarily")
    print(f"because the common cause failure rate alone ({lambda_ccf_fit:.3f} FIT) exceeds the target.")

if __name__ == '__main__':
    calculate_voter_failure_rate()
