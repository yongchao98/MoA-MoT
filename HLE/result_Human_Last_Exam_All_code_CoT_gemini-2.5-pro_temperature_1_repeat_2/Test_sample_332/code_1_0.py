import sys

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2oo3 system to meet ASIL C.
    """
    # --- Input Parameters ---
    lambda_ods1_fit = 500  # FIT
    lambda_ods2_fit = 400  # FIT
    lambda_ods3_fit = 700  # FIT
    beta = 0.10             # Common cause factor
    lambda_target_fit = 100 # ASIL C target in FIT
    t_mission_h = 10000     # Mission time in hours
    fit_to_per_hour = 1e-9  # Conversion factor from FIT to failures per hour

    print("Step-by-step calculation for the required voter failure rate.\n")
    print(f"The total system failure rate must be less than the target: λ_target = {lambda_target_fit} FIT.")
    print("The system failure rate is modeled as: λ_system = λ_2oo3_independent + λ_ccf + λ_voter\n")

    # --- 1. Calculate Common Cause Failure (CCF) Rate ---
    lambda_ods_avg_fit = (lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit) / 3.0
    lambda_ccf_fit = beta * lambda_ods_avg_fit
    
    # --- 2. Calculate Independent 2-out-of-3 Failure Rate ---
    # First, find the independent failure rates for each ODS
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit

    # Convert independent FIT rates to per-hour rates for the formula
    lambda_1_ind_h = lambda_ods1_ind_fit * fit_to_per_hour
    lambda_2_ind_h = lambda_ods2_ind_fit * fit_to_per_hour
    lambda_3_ind_h = lambda_ods3_ind_fit * fit_to_per_hour

    # Use the approximation formula for the average failure rate of a 2oo3 system over a mission time
    # λ_2oo3_avg ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * T_mission
    sum_of_products = (lambda_1_ind_h * lambda_2_ind_h +
                       lambda_1_ind_h * lambda_3_ind_h +
                       lambda_2_ind_h * lambda_3_ind_h)
    
    lambda_2oo3_ind_avg_h = sum_of_products * t_mission_h

    # Convert the result back to FIT
    lambda_2oo3_ind_avg_fit = lambda_2oo3_ind_avg_h / fit_to_per_hour

    # --- 3. Calculate the maximum allowed voter failure rate ---
    # λ_voter < λ_target - λ_ccf - λ_2oo3_independent
    lambda_voter_max_fit = lambda_target_fit - lambda_ccf_fit - lambda_2oo3_ind_avg_fit

    print("To find the maximum voter failure rate, we use the equation:")
    print("λvoter < λ_target - λ_ccf - λ_2oo3_independent\n")
    print("Calculated values for the equation:")
    print(f"  λ_target = {lambda_target_fit:.2f} FIT (The ASIL C requirement)")
    print(f"  λ_ccf = {lambda_ccf_fit:.2f} FIT (The contribution from common cause failures)")
    print(f"  λ_2oo3_independent = {lambda_2oo3_ind_avg_fit:.2f} FIT (The contribution from independent ODS failures)\n")
    
    print("Final result:")
    # The final print statement is formatted to match the user's request.
    print(f"λvoter < {lambda_target_fit:.2f} - {lambda_ccf_fit:.2f} - {lambda_2oo3_ind_avg_fit:.2f} FIT")
    print(f"λvoter < {lambda_voter_max_fit:.2f} FIT")

    # Capture the final numerical answer for the platform.
    # The output stream will not show this part.
    sys.stdout = open('/dev/null', 'w')
    print(f"<<<{lambda_voter_max_fit:.1f}>>>")
    sys.stdout = sys.__stdout__


if __name__ == "__main__":
    calculate_voter_failure_rate()
