def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet ASIL C requirements, considering common cause failures.
    """
    # Given failure rates in FIT (Failures In Time, 1 FIT = 1e-9 failures/hour)
    lambda_ods1_fit = 500
    lambda_ods2_fit = 400
    lambda_ods3_fit = 700

    # Common Cause Failure factor
    beta = 0.1

    # System requirement for ASIL C
    lambda_asil_c_fit = 100

    # Sample run time in hours
    t_hours = 10000

    # Conversion factor from FIT to failures per hour
    fit_to_per_hour = 1e-9

    # Step 1: Calculate the independent failure rate for each ODS
    # lambda_ind = (1 - beta) * lambda
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit

    # Convert independent failure rates from FIT to failures/hour for the DPF calculation
    l1_ind_h = lambda_ods1_ind_fit * fit_to_per_hour
    l2_ind_h = lambda_ods2_ind_fit * fit_to_per_hour
    l3_ind_h = lambda_ods3_ind_fit * fit_to_per_hour

    # Step 2: Calculate the Dual-Point Failure (DPF) rate of the ODS group
    # This is the rate of failure due to two independent ODS failures.
    # Formula: lambda_DPF = (l1*l2 + l1*l3 + l2*l3) * T
    lambda_dpf_h = (l1_ind_h * l2_ind_h + l1_ind_h * l3_ind_h + l2_ind_h * l3_ind_h) * t_hours
    
    # Convert the DPF rate back to FIT
    lambda_dpf_fit = lambda_dpf_h / fit_to_per_hour

    # Step 3: Calculate the Common Cause Failure (CCF) rate of the ODS group
    # Formula derived from beta = lambda_CCF / (lambda_DPF + lambda_CCF)
    lambda_ccf_fit = (beta * lambda_dpf_fit) / (1 - beta)

    # Step 4: Calculate the total failure rate of the ODS group
    lambda_ods_group_fit = lambda_dpf_fit + lambda_ccf_fit

    # Step 5: Calculate the maximum allowed failure rate for the voter
    # lambda_voter < lambda_ASIL_C - lambda_ODS_group
    lambda_voter_max_fit = lambda_asil_c_fit - lambda_ods_group_fit

    # Print the results and the final equation
    print("Calculation Steps:")
    print(f"1. Independent Failure Rates (FIT):")
    print(f"   λ_ODS1_ind = (1 - {beta}) * {lambda_ods1_fit} = {lambda_ods1_ind_fit:.3f} FIT")
    print(f"   λ_ODS2_ind = (1 - {beta}) * {lambda_ods2_fit} = {lambda_ods2_ind_fit:.3f} FIT")
    print(f"   λ_ODS3_ind = (1 - {beta}) * {lambda_ods3_fit} = {lambda_ods3_ind_fit:.3f} FIT")
    print("\n2. ODS Group Dual-Point Failure Rate (λ_DPF):")
    print(f"   λ_DPF = {lambda_dpf_fit:.3f} FIT")
    print("\n3. ODS Group Common Cause Failure Rate (λ_CCF):")
    print(f"   λ_CCF = ({beta} * {lambda_dpf_fit:.3f}) / (1 - {beta}) = {lambda_ccf_fit:.3f} FIT")
    print("\n4. Total ODS Group Failure Rate (λ_ODS_group):")
    print(f"   λ_ODS_group = λ_DPF + λ_CCF = {lambda_dpf_fit:.3f} + {lambda_ccf_fit:.3f} = {lambda_ods_group_fit:.3f} FIT")
    print("\n5. Required Voter Failure Rate (λ_voter):")
    print(f"   The system must be ASIL C compliant: λ_sys = λ_ODS_group + λ_voter < {lambda_asil_c_fit} FIT")
    print(f"   λ_voter < {lambda_asil_c_fit} - λ_ODS_group")
    print(f"   λ_voter < {lambda_asil_c_fit} - {lambda_ods_group_fit:.3f}")
    print(f"   λ_voter < {lambda_voter_max_fit:.3f} FIT")

    # Return the final numerical answer for the user
    return lambda_voter_max_fit

# Execute the calculation and print the result
final_answer = calculate_voter_failure_rate()
# The final answer is printed within the function.
# The problem asks for the answer in a specific format at the end.
# The value is rounded to 2 decimal places for the final answer format.
print(f"\n<<<λvoter < {final_answer:.2f} FIT>>>")