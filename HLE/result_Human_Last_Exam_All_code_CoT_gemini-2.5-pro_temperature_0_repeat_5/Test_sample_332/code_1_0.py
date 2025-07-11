def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    with common cause failures to meet a specific ASIL target.
    """
    # --- Given Data ---
    # Failure rates of the three ODS versions in FIT
    lambda_ods1 = 500.0
    lambda_ods2 = 400.0
    lambda_ods3 = 700.0

    # Common Cause Failure factor
    beta = 0.10

    # System target failure rate for ASIL C in FIT
    lambda_target = 100.0

    # Sample run time in hours
    T_mission = 10000.0

    # Conversion factor from FIT to failures per hour (1 FIT = 1e-9 failures/hour)
    FIT_TO_HOUR = 1e-9
    HOUR_TO_FIT = 1e9

    # --- Step 1: Calculate the Common Cause Failure (CCF) rate ---
    # For non-identical components, we use the average failure rate
    lambda_ods_avg = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3.0
    lambda_ccf = beta * lambda_ods_avg

    # --- Step 2: Calculate the failure rate of the 2oo3 system from independent failures ---
    # First, calculate the independent failure rate for each component
    lambda_ods1_ind = (1 - beta) * lambda_ods1
    lambda_ods2_ind = (1 - beta) * lambda_ods2
    lambda_ods3_ind = (1 - beta) * lambda_ods3

    # Convert independent FIT rates to per-hour rates for the formula
    l1_h = lambda_ods1_ind * FIT_TO_HOUR
    l2_h = lambda_ods2_ind * FIT_TO_HOUR
    l3_h = lambda_ods3_ind * FIT_TO_HOUR

    # Approximate the average failure rate for the 2oo3 redundant part
    # Formula: (λ1*λ2 + λ1*λ3 + λ2*λ3) * T
    lambda_2oo3_ind_h = (l1_h * l2_h + l1_h * l3_h + l2_h * l3_h) * T_mission

    # Convert the result back to FIT
    lambda_2oo3_ind = lambda_2oo3_ind_h * HOUR_TO_FIT

    # --- Step 3: Calculate the maximum allowed voter failure rate ---
    # The total failure rate is the sum of the parts: λ_target >= λ_2oo3_ind + λ_ccf + λ_voter
    lambda_voter_max = lambda_target - lambda_2oo3_ind - lambda_ccf

    # --- Step 4: Print the results ---
    print("The total system failure rate (λ_system) is the sum of the rates of its parts:")
    print("λ_system = λ_2oo3_independent + λ_common_cause + λ_voter")
    print("\nTo meet the ASIL C target, λ_system must be <= 100 FIT.")
    print("Therefore, the voter's failure rate must be:")
    print("λ_voter <= λ_target - λ_2oo3_independent - λ_common_cause")
    print("\nCalculated values:")
    print(f"λ_target = {lambda_target:.3f} FIT")
    print(f"λ_common_cause = {lambda_ccf:.3f} FIT")
    print(f"λ_2oo3_independent = {lambda_2oo3_ind:.3f} FIT")
    
    print("\nFinal Equation:")
    # Using the format "λvoter < X FIT" as requested
    print(f"λ_voter < {lambda_target:.3f} - {lambda_2oo3_ind:.3f} - {lambda_ccf:.3f} FIT")
    print(f"λ_voter < {lambda_voter_max:.3f} FIT")
    
    # Return the final numerical answer for the wrapper
    return lambda_voter_max

# Execute the calculation and print the result
final_answer = calculate_voter_failure_rate()
print(f"\n<<< {final_answer:.3f} >>>")
