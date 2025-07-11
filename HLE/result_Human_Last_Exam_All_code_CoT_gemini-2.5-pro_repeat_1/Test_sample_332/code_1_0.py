import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate of the voter for the system to be ASIL C compliant.
    """
    # --- Step 0: Define Initial Parameters ---
    lambda_ods1_fit = 500  # FIT
    lambda_ods2_fit = 400  # FIT
    lambda_ods3_fit = 700  # FIT
    beta = 0.10             # Common cause factor
    t_mission = 10000       # Sample run time in hours
    lambda_system_target_fit = 100 # ASIL C target in FIT
    
    # FIT is failures per 10^9 hours
    fit_to_rate_per_hour = 1e-9

    print("--- System Parameters ---")
    print(f"ODS1 Failure Rate (λ_ODS1): {lambda_ods1_fit} FIT")
    print(f"ODS2 Failure Rate (λ_ODS2): {lambda_ods2_fit} FIT")
    print(f"ODS3 Failure Rate (λ_ODS3): {lambda_ods3_fit} FIT")
    print(f"Common Cause Factor (β): {beta}")
    print(f"Mission Time (t): {t_mission} hours")
    print(f"System Target Failure Rate (λ_target): {lambda_system_target_fit} FIT\n")

    # --- Step 1: Calculate Independent Failure Rates ---
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit

    print("--- Step 1: Independent Failure Rate Calculation ---")
    print(f"Independent λ_ODS1 = (1 - {beta}) * {lambda_ods1_fit} = {lambda_ods1_ind_fit} FIT")
    print(f"Independent λ_ODS2 = (1 - {beta}) * {lambda_ods2_fit} = {lambda_ods2_ind_fit} FIT")
    print(f"Independent λ_ODS3 = (1 - {beta}) * {lambda_ods3_fit} = {lambda_ods3_ind_fit} FIT\n")

    # Convert independent FIT rates to per-hour rates for the formula
    l1_ind_hr = lambda_ods1_ind_fit * fit_to_rate_per_hour
    l2_ind_hr = lambda_ods2_ind_fit * fit_to_rate_per_hour
    l3_ind_hr = lambda_ods3_ind_fit * fit_to_rate_per_hour

    # --- Step 2: Calculate 2-out-of-3 Independent System Failure Rate ---
    # λ_2oo3 ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * t
    lambda_2oo3_ind_hr = (l1_ind_hr * l2_ind_hr + 
                           l1_ind_hr * l3_ind_hr + 
                           l2_ind_hr * l3_ind_hr) * t_mission
    
    # Convert back to FIT
    lambda_2oo3_ind_fit = lambda_2oo3_ind_hr / fit_to_rate_per_hour
    
    print("--- Step 2: 2-out-of-3 ODS Group Failure Rate (Independent) ---")
    print(f"λ_independent_2oo3 ≈ (λ_ind1*λ_ind2 + λ_ind1*λ_ind3 + λ_ind2*λ_ind3) * t")
    print(f"λ_independent_2oo3 = {lambda_2oo3_ind_fit:.4f} FIT\n")

    # --- Step 3: Calculate Common Cause Failure Rate ---
    lambda_avg_fit = (lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit) / 3
    lambda_ccf_fit = beta * lambda_avg_fit

    print("--- Step 3: Common Cause Failure (CCF) Rate Calculation ---")
    print(f"Average ODS Failure Rate (λ_avg) = ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) / 3 = {lambda_avg_fit:.4f} FIT")
    print(f"CCF Rate (λ_ccf) = β * λ_avg = {beta} * {lambda_avg_fit:.4f} = {lambda_ccf_fit:.4f} FIT\n")

    # --- Step 4: Calculate Required Voter Failure Rate ---
    # λ_voter < λ_target - (λ_independent_2oo3 + λ_ccf)
    lambda_ods_group_fit = lambda_2oo3_ind_fit + lambda_ccf_fit
    lambda_voter_max_fit = lambda_system_target_fit - lambda_ods_group_fit

    print("--- Step 4: Final Voter Failure Rate Calculation ---")
    print(f"Total ODS Group Failure Rate = λ_independent_2oo3 + λ_ccf = {lambda_2oo3_ind_fit:.4f} + {lambda_ccf_fit:.4f} = {lambda_ods_group_fit:.4f} FIT")
    print(f"Required λ_voter < λ_target - λ_ODS_group")
    print(f"Required λ_voter < {lambda_system_target_fit} - ({lambda_2oo3_ind_fit:.4f} + {lambda_ccf_fit:.4f})")
    print(f"Required λ_voter < {lambda_system_target_fit} - {lambda_ods_group_fit:.4f}")
    print(f"Required λ_voter < {lambda_voter_max_fit:.4f} FIT\n")
    
    # Final Answer
    print(f"Final Answer: The required failure rate for the voter is:")
    print(f"λvoter < {lambda_voter_max_fit:.4f} FIT")
    
    # Returning the final number for the user's format.
    return lambda_voter_max_fit

# Execute the calculation and print the result
final_answer = calculate_voter_failure_rate()
# The final answer format as requested by the user prompt
print(f'<<<{final_answer:.4f}>>>')
