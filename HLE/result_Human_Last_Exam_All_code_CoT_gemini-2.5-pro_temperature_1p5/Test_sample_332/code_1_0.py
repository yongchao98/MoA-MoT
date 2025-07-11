import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet a specific ASIL target, considering common cause failures.
    """
    # Given values from the problem description
    lambda_ods1_fit = 500
    lambda_ods2_fit = 400
    lambda_ods3_fit = 700
    beta = 0.10
    lambda_system_target_fit = 100
    t_hours = 10000  # Mission time in hours
    fit_to_per_hour = 1e-9

    print("System Reliability Calculation")
    print("---------------------------------")
    print("The total system failure rate is the sum of the components in series:")
    print("λ_system = λ_2oo3_independent + λ_ccf + λ_voter\n")
    print(f"The target system failure rate for ASIL C compliance is λ_target = {lambda_system_target_fit} FIT.")
    
    # --- Step 1: Calculate Common Cause Failure (CCF) rate (λ_ccf) ---
    # This is modeled as β times the average failure rate of the ODS components.
    lambda_sum_fit = lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit
    lambda_avg_fit = lambda_sum_fit / 3
    lambda_ccf_fit = beta * lambda_avg_fit
    
    print("\nStep 1: Calculate the Common Cause Failure (CCF) rate (λ_ccf)")
    print(f"λ_ccf = β * (λ_ODS1 + λ_ODS2 + λ_ODS3) / 3")
    print(f"λ_ccf = {beta} * ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) / 3 = {lambda_ccf_fit:.2f} FIT")

    # --- Step 2: Calculate the failure rate of the 2oo3 system from independent failures (λ_2oo3_independent) ---
    # First, find the independent failure rates for each ODS
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit
    
    # Convert independent FIT rates to per-hour rates
    lambda_1_ind_hr = lambda_ods1_ind_fit * fit_to_per_hour
    lambda_2_ind_hr = lambda_ods2_ind_fit * fit_to_per_hour
    lambda_3_ind_hr = lambda_ods3_ind_fit * fit_to_per_hour
    
    # Calculate the probability of failure for each component over mission time t
    F1_ind = 1 - math.exp(-lambda_1_ind_hr * t_hours)
    F2_ind = 1 - math.exp(-lambda_2_ind_hr * t_hours)
    F3_ind = 1 - math.exp(-lambda_3_ind_hr * t_hours)
    
    # Calculate the probability of failure for the 2oo3 redundant system (due to independent failures)
    # The system fails if at least two components fail.
    # P(sys_fail) = P(F1∩F2) + P(F1∩F3) + P(F2∩F3) - 2*P(F1∩F2∩F3)
    F_2oo3_ind = (F1_ind * F2_ind) + (F1_ind * F3_ind) + (F2_ind * F3_ind) - 2 * (F1_ind * F2_ind * F3_ind)

    # Calculate the average failure rate (in FIT) over the mission time t
    lambda_2oo3_ind_avg_hr = F_2oo3_ind / t_hours
    lambda_2oo3_ind_avg_fit = lambda_2oo3_ind_avg_hr / fit_to_per_hour
    
    print("\nStep 2: Calculate the equivalent failure rate from independent failures (λ_2oo3_independent)")
    print(f"Based on a mission time of {t_hours} hours, the average rate of the 2oo3 subsystem is calculated.")
    print(f"λ_2oo3_independent = {lambda_2oo3_ind_avg_fit:.2f} FIT")

    # --- Step 3: Calculate the required voter failure rate (λ_voter) ---
    # From λ_system = λ_2oo3_independent + λ_ccf + λ_voter
    # λ_voter = λ_system_target - λ_ccf - λ_2oo3_independent
    lambda_voter_req_fit = lambda_system_target_fit - lambda_ccf_fit - lambda_2oo3_ind_avg_fit
    
    print("\nStep 3: Calculate the maximum allowed failure rate for the voter (λ_voter)")
    print("λ_voter < λ_target - λ_ccf - λ_2oo3_independent")
    # Outputting each number in the final equation
    print(f"λ_voter < {lambda_system_target_fit} - {lambda_ccf_fit:.2f} - {lambda_2oo3_ind_avg_fit:.2f}")
    print(f"λ_voter < {lambda_voter_req_fit:.2f} FIT")
    
    return lambda_voter_req_fit

if __name__ == '__main__':
    result = calculate_voter_failure_rate()
    # The final answer format as requested
    # print(f"\n<<<λvoter < {result:.2f} FIT>>>")
    # print(f"<<<{result:.2f}>>>")