import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet a specific ASIL target, considering common cause failures.
    """
    # --- Inputs from the problem statement ---
    lambda_ods1_fit = 500.0
    lambda_ods2_fit = 400.0
    lambda_ods3_fit = 700.0
    beta = 0.10
    m = 3  # Number of ODS units in the common cause group
    t_mission_h = 10000.0
    lambda_system_target_fit = 100.0

    # --- Unit conversion factors ---
    fit_to_per_hour = 1e-9
    per_hour_to_fit = 1e9

    print("Step-by-step calculation for the required voter failure rate (λ_voter).\n")
    print(f"System Target: The total system failure rate must be < {lambda_system_target_fit} FIT (ASIL C).")
    print("System Model: λ_system = λ_voter + λ_common_cause + λ_dual_point < 100 FIT\n")

    # 1. Calculate the Common Cause Failure (CCF) rate
    sum_lambda_ods_fit = lambda_ods1_fit + lambda_ods2_fit + lambda_ods3_fit
    # For a group of m components, the rate of a single CCF event is modeled as (β/m) * Σλ
    lambda_ccf_fit = (beta / m) * sum_lambda_ods_fit
    
    print("1. Calculate the Common Cause Failure rate (λ_ccf):")
    print(f"   λ_ccf = (β / m) * (λ_ODS1 + λ_ODS2 + λ_ODS3)")
    print(f"   λ_ccf = ({beta} / {m}) * ({lambda_ods1_fit} + {lambda_ods2_fit} + {lambda_ods3_fit}) FIT")
    print(f"   λ_ccf = {lambda_ccf_fit:.2f} FIT\n")

    # 2. Calculate the Dual-Point Failure (DPF) rate from independent failures
    lambda_ods1_ind_fit = (1 - beta) * lambda_ods1_fit
    lambda_ods2_ind_fit = (1 - beta) * lambda_ods2_fit
    lambda_ods3_ind_fit = (1 - beta) * lambda_ods3_fit

    # Convert independent FIT rates to per-hour for the formula
    l1_ind_ph = lambda_ods1_ind_fit * fit_to_per_hour
    l2_ind_ph = lambda_ods2_ind_fit * fit_to_per_hour
    l3_ind_ph = lambda_ods3_ind_fit * fit_to_per_hour

    # The average failure rate for a 2oo3 system is approximated by:
    # λ_dpf ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * t_mission
    sum_of_products_ph2 = (l1_ind_ph * l2_ind_ph) + \
                          (l1_ind_ph * l3_ind_ph) + \
                          (l2_ind_ph * l3_ind_ph)
    lambda_dpf_ph = sum_of_products_ph2 * t_mission_h
    lambda_dpf_fit = lambda_dpf_ph * per_hour_to_fit
    
    print("2. Calculate the Dual-Point Failure rate (λ_dpf) for the 2-out-of-3 architecture:")
    print(f"   Independent failure rates (FIT): λ_ind1={lambda_ods1_ind_fit:.0f}, λ_ind2={lambda_ods2_ind_fit:.0f}, λ_ind3={lambda_ods3_ind_fit:.0f}")
    print(f"   Using the approximation for a mission time of {t_mission_h:.0f}h:")
    print(f"   λ_dpf = {lambda_dpf_fit:.2f} FIT\n")

    # 3. Calculate the available budget for the voter
    lambda_ods_subsystem_fit = lambda_ccf_fit + lambda_dpf_fit
    max_lambda_voter_fit = lambda_system_target_fit - lambda_ods_subsystem_fit
    
    print("3. Calculate the maximum allowed failure rate for the voter (λ_voter):")
    print(f"   λ_voter < λ_system_target - (λ_ccf + λ_dpf)")
    print(f"   λvoter < {lambda_system_target_fit:.2f} - ({lambda_ccf_fit:.2f} + {lambda_dpf_fit:.2f})")
    print(f"   λvoter < {lambda_system_target_fit:.2f} - {lambda_ods_subsystem_fit:.2f}")
    print("\nFinal Result:")
    print(f"λvoter < {max_lambda_voter_fit:.2f} FIT")
    
    # Final answer in the required format
    print(f"\n<<< {max_lambda_voter_fit:.2f} >>>")

if __name__ == '__main__':
    calculate_voter_failure_rate()