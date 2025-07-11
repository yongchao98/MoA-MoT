import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate of a voter in a 2-out-of-3 system
    to meet a specific ASIL target.
    """
    # --- Given Parameters ---
    # Failure rates of the three ODS units in FIT (Failures In Time, 10^-9 failures/hour)
    lambda_ods1 = 500
    lambda_ods2 = 400
    lambda_ods3 = 700

    # Common Cause Failure (CCF) Beta factor
    beta = 0.10

    # ASIL C target failure rate for the complete system in FIT
    lambda_target = 100

    # Sample run time in hours for reliability calculation
    mission_time_h = 10000

    # --- Step 1: Calculate the Common Cause Failure (CCF) Rate ---
    # For non-identical components, a common approach is to use the average failure rate.
    ods_failure_rates = [lambda_ods1, lambda_ods2, lambda_ods3]
    avg_lambda = sum(ods_failure_rates) / len(ods_failure_rates)
    lambda_ccf = beta * avg_lambda

    # --- Step 2: Calculate the Independent Failure Rate for the 2oo3 System ---
    # The independent portion of the failure rate for each ODS
    lambda_ods1_ind = (1 - beta) * lambda_ods1
    lambda_ods2_ind = (1 - beta) * lambda_ods2
    lambda_ods3_ind = (1 - beta) * lambda_ods3

    # Conversion factor from FIT to failures per hour
    fit_to_per_hour = 1e-9

    # Calculate the failure probability F(t) = 1 - R(t) ≈ λ*t for each component over the mission time
    F1_ind = lambda_ods1_ind * fit_to_per_hour * mission_time_h
    F2_ind = lambda_ods2_ind * fit_to_per_hour * mission_time_h
    F3_ind = lambda_ods3_ind * fit_to_per_hour * mission_time_h

    # The unreliability (failure probability) of a 2oo3 system is: F1*F2 + F1*F3 + F2*F3 - 2*F1*F2*F3
    F_2oo3_ind = (F1_ind * F2_ind) + (F1_ind * F3_ind) + (F2_ind * F3_ind) - (2 * F1_ind * F2_ind * F3_ind)

    # Convert the system's failure probability back to an average failure rate in FIT
    lambda_independent_2oo3_per_hour = F_2oo3_ind / mission_time_h
    lambda_independent_2oo3 = lambda_independent_2oo3_per_hour / fit_to_per_hour

    # --- Step 3: Calculate the Maximum Permissible Voter Failure Rate ---
    # The total system failure rate is the sum of the components in series:
    # λ_target = λ_independent_2oo3 + λ_ccf + λ_voter
    # Therefore, λ_voter = λ_target - λ_independent_2oo3 - λ_ccf
    lambda_voter_max = lambda_target - lambda_independent_2oo3 - lambda_ccf

    # --- Step 4: Output the Final Result ---
    print("The total system failure rate (λ_system) must be <= λ_target.")
    print("λ_system = λ_independent_2oo3 + λ_ccf + λ_voter\n")
    print("Calculation steps:")
    print(f"1. Common Cause Failure Rate (λ_ccf):")
    print(f"   λ_ccf = β * λ_average = {beta} * (({lambda_ods1} + {lambda_ods2} + {lambda_ods3}) / 3) = {lambda_ccf:.2f} FIT")
    print(f"2. Independent 2-out-of-3 Failure Rate (λ_independent_2oo3): {lambda_independent_2oo3:.2f} FIT")
    print("\nDeriving the required voter failure rate:")
    print(f"λ_voter <= λ_target - λ_independent_2oo3 - λ_ccf")
    print(f"λ_voter <= {lambda_target} FIT - {lambda_independent_2oo3:.2f} FIT - {lambda_ccf:.2f} FIT")
    print(f"λ_voter <= {lambda_voter_max:.2f} FIT")


# Execute the calculation
calculate_voter_failure_rate()