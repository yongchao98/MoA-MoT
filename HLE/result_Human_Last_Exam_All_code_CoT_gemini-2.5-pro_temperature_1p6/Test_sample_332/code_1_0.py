def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    with common cause failures to meet a specific ASIL target.
    """
    # --- Step 1: Define Constants ---
    lambda_ods1_fit = 500  # FIT
    lambda_ods2_fit = 400  # FIT
    lambda_ods3_fit = 700  # FIT
    beta = 0.1
    t_mission_h = 10000  # hours
    lambda_system_target_fit = 100  # ASIL C Target in FIT
    
    fit_to_per_hour = 1e-9

    # --- Step 2: Calculate Independent Failure Rates ---
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit

    # --- Step 3: Calculate the 2oo3 Independent Failure Rate ---
    # Convert independent FIT rates to per-hour rates for the formula
    l1 = lambda_ods1_ind_fit * fit_to_per_hour
    l2 = lambda_ods2_ind_fit * fit_to_per_hour
    l3 = lambda_ods3_ind_fit * fit_to_per_hour
    
    # Formula for average PFH of a 2oo3 system with non-identical components
    # PFH_avg = T_mission * (λ1λ2 + λ1λ3 + λ2λ3)
    lambda_2oo3_ind_per_hour = t_mission_h * (l1 * l2 + l1 * l3 + l2 * l3)
    
    # Convert the result back to FIT
    lambda_2oo3_ind_fit = lambda_2oo3_ind_per_hour / fit_to_per_hour

    # --- Step 4: Calculate the Common Cause Failure Rate ---
    # We assume beta is applied to the average failure rate of the ODS units.
    lambda_ods_avg_fit = (lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit) / 3
    lambda_ccf_fit = beta * lambda_ods_avg_fit
    
    # --- Step 5: Calculate the Maximum Voter Failure Rate ---
    # λ_target > λ_2oo3_ind + λ_ccf + λ_voter
    # λ_voter < λ_target - λ_2oo3_ind - λ_ccf
    lambda_voter_max_fit = lambda_system_target_fit - lambda_2oo3_ind_fit - lambda_ccf_fit
    
    # --- Final Output ---
    print("The final calculation is based on the equation:")
    print("λ_voter < λ_target - λ_2oo3_independent - λ_common_cause")
    print("\nSubstituting the calculated values:")
    # Using format specifiers to control the number of decimal places for clarity
    print("λ_voter < {:.3f} FIT - {:.3f} FIT - {:.3f} FIT".format(
        lambda_system_target_fit,
        lambda_2oo3_ind_fit,
        lambda_ccf_fit
    ))
    
    print("\nResult:")
    print("λvoter < {:.3f} FIT".format(lambda_voter_max_fit))


calculate_voter_failure_rate()
<<<39.944>>>