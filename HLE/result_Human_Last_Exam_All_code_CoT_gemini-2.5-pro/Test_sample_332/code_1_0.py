import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    with common cause failures to meet a specific ASIL target.
    """
    # Given parameters from the problem description
    lambda_ods1 = 500  # FIT
    lambda_ods2 = 400  # FIT
    lambda_ods3 = 700  # FIT
    beta = 0.1
    lambda_target = 100  # ASIL C target in FIT
    t_mission = 10000  # Mission time in hours

    print("This script calculates the required voter failure rate for a 2-out-of-3 system.")
    print("-" * 60)
    print("Step 1: Define the total system failure rate equation.")
    print(f"The system must meet the target: λ_system <= {lambda_target} FIT")
    print("The total failure rate is a sum of its parts:")
    print("λ_system = λ_mpf_independent + λ_ccf + λ_voter")
    print("-" * 60)

    # Step 2: Calculate the Common Cause Failure (CCF) rate (λ_ccf).
    # For non-identical components, we use the average failure rate.
    lambda_avg = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
    lambda_ccf = beta * lambda_avg

    print("Step 2: Calculate the Common Cause Failure (CCF) rate (λ_ccf).")
    print(f"   λ_avg = ({lambda_ods1} + {lambda_ods2} + {lambda_ods3}) / 3 = {lambda_avg:.3f} FIT")
    print(f"   λ_ccf = β * λ_avg = {beta} * {lambda_avg:.3f} = {lambda_ccf:.3f} FIT")
    print("-" * 60)

    # Step 3: Calculate the Multiple-Point Failure rate from independent failures (λ_mpf_ind).
    # First, get the independent failure rate for each ODS.
    lambda_ods1_ind = lambda_ods1 * (1 - beta)
    lambda_ods2_ind = lambda_ods2 * (1 - beta)
    lambda_ods3_ind = lambda_ods3 * (1 - beta)

    # Convert FIT (10^-9 failures/hr) to failures/hr for the formula.
    fit_to_per_hour = 1e-9
    lambda_ods1_ind_h = lambda_ods1_ind * fit_to_per_hour
    lambda_ods2_ind_h = lambda_ods2_ind * fit_to_per_hour
    lambda_ods3_ind_h = lambda_ods3_ind * fit_to_per_hour

    # Approximate the average failure rate for dual-point faults over the mission time.
    lambda_mpf_ind_h = (lambda_ods1_ind_h * lambda_ods2_ind_h +
                       lambda_ods1_ind_h * lambda_ods3_ind_h +
                       lambda_ods2_ind_h * lambda_ods3_ind_h) * t_mission

    # Convert the result back from failures/hr to FIT.
    lambda_mpf_ind = lambda_mpf_ind_h / fit_to_per_hour

    print("Step 3: Calculate the Multiple-Point Failure rate (λ_mpf_independent).")
    print("   First, find the independent failure rate of each component:")
    print(f"   λ_ODS1_ind = {lambda_ods1} * (1 - {beta}) = {lambda_ods1_ind:.1f} FIT")
    print(f"   λ_ODS2_ind = {lambda_ods2} * (1 - {beta}) = {lambda_ods2_ind:.1f} FIT")
    print(f"   λ_ODS3_ind = {lambda_ods3} * (1 - {beta}) = {lambda_ods3_ind:.1f} FIT")
    print("\n   Then, calculate the system's multiple-point failure rate from these:")
    print(f"   λ_mpf_ind ≈ (λ1_ind*λ2_ind + λ1_ind*λ3_ind + λ2_ind*λ3_ind) * T_mission")
    print(f"   λ_mpf_ind = {lambda_mpf_ind:.3f} FIT")
    print("-" * 60)

    # Step 4: Calculate the maximum allowed failure rate for the voter.
    lambda_voter_max = lambda_target - lambda_mpf_ind - lambda_ccf

    print("Step 4: Calculate the required voter failure rate (λ_voter).")
    print("   λ_voter < λ_target - λ_mpf_independent - λ_ccf")
    print(f"   λ_voter < {lambda_target} - {lambda_mpf_ind:.3f} - {lambda_ccf:.3f}")
    final_result_str = f"λ_voter < {lambda_voter_max:.3f} FIT"
    print(f"   {final_result_str}")
    print("-" * 60)

if __name__ == '__main__':
    calculate_voter_failure_rate()