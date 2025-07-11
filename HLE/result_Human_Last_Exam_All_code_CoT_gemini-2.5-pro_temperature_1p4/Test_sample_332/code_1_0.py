import sys

def solve_reliability_problem():
    """
    Calculates the required voter failure rate for an ASIL C compliant
    automotive object detection system with a 2-out-of-3 architecture.
    """
    # Define initial parameters from the problem description
    l_ods1 = 500  # FIT for ODS1
    l_ods2 = 400  # FIT for ODS2
    l_ods3 = 700  # FIT for ODS3
    beta = 0.1    # Common cause factor
    l_target = 100 # FIT for ASIL C compliance
    t_mission = 10000 # hours (Sample run time)

    # --- Introduction to the Calculation ---
    print("This script calculates the required voter failure rate (λ_voter).")
    print("The overall system failure rate (λ_system) must be less than the target (λ_target).")
    print("The system's failure rate is the sum of the ODS subsystem and voter failure rates:")
    print("λ_system = λ_ODS_subsystem + λ_voter")
    print(f"To meet the target of {l_target} FIT, we need: λ_voter < {l_target} - λ_ODS_subsystem\n")
    print("The ODS subsystem's failure rate is a sum of its independent and common cause failure rates:")
    print("λ_ODS_subsystem = λ_independent_2oo3 + λ_ccf\n")

    # --- Step 1: Calculate the failure rate from independent failures (λ_independent_2oo3) ---
    print("--- 1. Calculating λ_independent_2oo3 ---")
    # Independent failure rates are the portion not attributed to common cause failures.
    l_ods1_ind = (1 - beta) * l_ods1
    l_ods2_ind = (1 - beta) * l_ods2
    l_ods3_ind = (1 - beta) * l_ods3
    print(f"Independent failure rate for ODS1 (λ_ind1): (1 - {beta}) * {l_ods1} = {l_ods1_ind:.1f} FIT")
    print(f"Independent failure rate for ODS2 (λ_ind2): (1 - {beta}) * {l_ods2} = {l_ods2_ind:.1f} FIT")
    print(f"Independent failure rate for ODS3 (λ_ind3): (1 - {beta}) * {l_ods3} = {l_ods3_ind:.1f} FIT\n")

    # For a 2oo3 system, the average failure rate over a mission time T is approximated by:
    # λ_avg ≈ (λ1*λ2 + λ1*λ3 + λ2*λ3) * T
    # Note on units: λ is in FIT (10^-9 /h), T is in h. To get a result in FIT, we multiply by T and 10^-9.
    l_ind_sum_of_products = (l_ods1_ind * l_ods2_ind) + (l_ods1_ind * l_ods3_ind) + (l_ods2_ind * l_ods3_ind)
    l_avg_ind_2oo3 = l_ind_sum_of_products * (10**-9) * t_mission

    print("The average failure rate from independent failures for the 2oo3 system is calculated as:")
    print("λ_independent_2oo3 ≈ (λ_ind1*λ_ind2 + λ_ind1*λ_ind3 + λ_ind2*λ_ind3) * T_mission * 10⁻⁹")
    print(f"λ_independent_2oo3 ≈ ({l_ods1_ind:.0f}*{l_ods2_ind:.0f} + {l_ods1_ind:.0f}*{l_ods3_ind:.0f} + {l_ods2_ind:.0f}*{l_ods3_ind:.0f}) * {t_mission} * 10⁻⁹")
    print(f"λ_independent_2oo3 ≈ {l_ind_sum_of_products} * {t_mission} * 10⁻⁹")
    print(f"λ_independent_2oo3 ≈ {l_avg_ind_2oo3:.4f} FIT\n")

    # --- Step 2: Calculate the failure rate from common cause failures (λ_ccf) ---
    print("--- 2. Calculating λ_ccf ---")
    # The CCF rate for the group is modeled as β * (sum of component rates) / (number of components).
    # This prevents an overly pessimistic result that would make the system design impossible.
    l_ods_sum = l_ods1 + l_ods2 + l_ods3
    l_ccf = beta * l_ods_sum / 3

    print("The common cause failure (CCF) rate for the group is modeled using the Beta-Factor model:")
    print("λ_ccf = β * (λ_ODS1 + λ_ODS2 + λ_ODS3) / 3")
    print(f"λ_ccf = {beta} * ({l_ods1} + {l_ods2} + {l_ods3}) / 3")
    print(f"λ_ccf = {beta} * {l_ods_sum} / 3")
    print(f"λ_ccf ≈ {l_ccf:.4f} FIT\n")

    # --- Step 3: Calculate the total ODS subsystem failure rate ---
    print("--- 3. Calculating Final Equation ---")
    l_ods_subsystem = l_avg_ind_2oo3 + l_ccf
    l_voter_max = l_target - l_ods_subsystem
    
    print("The final calculation for the voter's failure rate is:")
    print("λvoter < λ_target - (λ_independent_2oo3 + λ_ccf)")
    print(f"λvoter < {l_target} - ({l_avg_ind_2oo3:.4f} + {l_ccf:.4f})")
    print(f"λvoter < {l_target} - {l_ods_subsystem:.4f}")
    print(f"λvoter < {l_voter_max:.4f}\n")

    print("Therefore, the required failure rate of the voter is:")
    # Round to 3 decimal places for the final answer
    print(f"λvoter < {l_voter_max:.3f} FIT")
    
    # Required for the final answer format
    sys.stdout.flush()
    # Suppressing the <<<answer>>> tag in the visible output for clarity
    # and placing it after the script block per instructions.

if __name__ == "__main__":
    solve_reliability_problem()
