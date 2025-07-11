import math

def solve_voter_failure_rate():
    """
    Calculates the required failure rate of a voter in a 2-out-of-3 system
    with common cause failures to meet a specific ASIL target.
    """
    # Step 1: Define system parameters
    lambda_ods1 = 500  # FIT
    lambda_ods2 = 400  # FIT
    lambda_ods3 = 700  # FIT
    beta = 0.10  # Common Cause Factor
    lambda_system_target = 100  # FIT for ASIL C
    t_mission = 10000  # hours

    print("Step 1: System Parameters")
    print(f"λODS1 = {lambda_ods1} FIT")
    print(f"λODS2 = {lambda_ods2} FIT")
    print(f"λODS3 = {lambda_ods3} FIT")
    print(f"Common Cause Factor (β) = {beta}")
    print(f"System Target Failure Rate (λ_target) = {lambda_system_target} FIT")
    print(f"Mission Time (t) = {t_mission} h\n")

    # Step 2: Calculate the failure rate of the ODS subsystem (λ_ODS_sys)

    # a) Calculate independent failure rates
    lambda_ods1_ind = (1 - beta) * lambda_ods1
    lambda_ods2_ind = (1 - beta) * lambda_ods2
    lambda_ods3_ind = (1 - beta) * lambda_ods3

    # b) Calculate the failure rate from independent failures in the 2oo3 configuration
    # Formula: λ_ind_2oo3 ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * t * 10⁻⁹
    products_sum = (lambda_ods1_ind * lambda_ods2_ind +
                    lambda_ods1_ind * lambda_ods3_ind +
                    lambda_ods2_ind * lambda_ods3_ind)
    lambda_ind_2oo3 = products_sum * t_mission * 1e-9

    # c) Calculate the failure rate from Common Cause Failures (CCF)
    # Since components are not identical, we use the average failure rate
    lambda_ods_avg = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
    lambda_ccf = beta * lambda_ods_avg

    # d) Calculate the total failure rate of the ODS subsystem
    lambda_ods_system = lambda_ind_2oo3 + lambda_ccf

    print("Step 2: Calculate ODS Subsystem Failure Rate (λ_ODS_sys)")
    print(f"   a) The failure rate from independent failures in the 2oo3 configuration (λ_ind_2oo3) is {lambda_ind_2oo3:.3f} FIT.")
    print(f"   b) The failure rate from common cause failures (λ_ccf) is {lambda_ccf:.3f} FIT.")
    print(f"   c) The total ODS subsystem failure rate is λ_ODS_sys = λ_ind_2oo3 + λ_ccf")
    print(f"      λ_ODS_sys = {lambda_ind_2oo3:.3f} + {lambda_ccf:.3f} = {lambda_ods_system:.3f} FIT\n")

    # Step 3: Calculate the required voter failure rate (λ_voter)
    # The total system failure rate is λ_sys = λ_ODS_sys + λ_voter
    # To meet the target: λ_ODS_sys + λ_voter < λ_target
    lambda_voter_max = lambda_system_target - lambda_ods_system

    print("Step 3: Calculate the required voter failure rate (λ_voter)")
    print("The final equation is: λ_voter < λ_target - λ_ODS_sys")
    print(f"λ_voter < {lambda_system_target} - {lambda_ods_system:.3f}")
    print(f"λ_voter < {lambda_voter_max:.3f} FIT\n")
    
    print("Therefore, the required failure rate of the voter must be less than this value.")
    
    # Return the final numeric answer for the <<<>>> format
    return lambda_voter_max

# Execute the function and capture the result
result = solve_voter_failure_rate()
# The final answer format requires the raw number to be enclosed in <<<>>>
# print(f"\n<<< {result:.3f} >>>") #This line is commented out as per instructions to not include multiple code blocks. The value will be manually placed.