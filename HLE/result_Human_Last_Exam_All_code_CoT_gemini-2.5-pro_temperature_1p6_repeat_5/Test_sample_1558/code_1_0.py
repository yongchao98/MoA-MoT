def analyze_hematopoiesis_data():
    """
    Analyzes and prints key findings from the described hematopoiesis experiments.
    """
    print("--- Data Analysis ---")

    # Experiment 1: RTI effect on RBC in pregnant mice
    print("\n[Experiment 1: Effect of RTI on Red Blood Cells in Pregnant Mice]")
    preg_control_rbc_exp1 = 10  # in millions per ul
    preg_rti_rbc_exp1 = 8       # in millions per ul
    print(f"The number of Red Blood Cells in pregnant control mice is {preg_control_rbc_exp1}x10^6 per ul.")
    print(f"The number of Red Blood Cells in pregnant RTI-treated mice is {preg_rti_rbc_exp1}x10^6 per ul.")
    print(f"Comparison Equation: {preg_control_rbc_exp1}x10^6 > {preg_rti_rbc_exp1}x10^6")
    print("Conclusion: Inhibiting transposable elements with RTI decreases the RBC count, suggesting that active transposable elements increase erythropoiesis during pregnancy.")

    # Experiment 2: STING deletion effect on RBC in pregnant mice
    print("\n[Experiment 2: Effect of STING Deletion on Red Blood Cells in Pregnant Mice]")
    preg_control_rbc_exp2 = 13  # in millions per ul
    preg_dsting_rbc_exp2 = 8    # in millions per ul
    print(f"The number of Red Blood Cells in pregnant control mice is {preg_control_rbc_exp2}x10^6 per ul.")
    print(f"The number of Red Blood Cells in pregnant delta STING mice is {preg_dsting_rbc_exp2}x10^6 per ul.")
    print(f"Comparison Equation: {preg_control_rbc_exp2}x10^6 > {preg_dsting_rbc_exp2}x10^6")
    print("Conclusion: Deleting the STING gene (part of the immune system) decreases the RBC count, suggesting that the STING-mediated immune response influences erythropoiesis.")

    # Experiment 3: IFNAR1 deletion effect on spleen progenitors in pregnant mice
    print("\n[Experiment 3: Effect of IFNAR1 Deletion on Spleen Progenitors in Pregnant Mice]")
    # HSC analysis
    preg_control_hsc = 0.003  # as percentage
    preg_difnar1_hsc = 0.002  # as percentage
    print(f"The percentage of HSC in spleen cells of pregnant control mice is {preg_control_hsc}%.")
    print(f"The percentage of HSC in spleen cells of pregnant delta IFNAR1 mice is {preg_difnar1_hsc}%.")
    print(f"Comparison Equation: {preg_control_hsc}% > {preg_difnar1_hsc}%")
    # MPP analysis
    preg_control_mpp = 0.004  # as percentage
    preg_difnar1_mpp = 0.002  # as percentage
    print(f"The percentage of MPP in spleen cells of pregnant control mice is {preg_control_mpp}%.")
    print(f"The percentage of MPP in spleen cells of pregnant delta IFNAR1 mice is {preg_difnar1_mpp}%.")
    print(f"Comparison Equation: {preg_control_mpp}% > {preg_difnar1_mpp}%")
    print("Conclusion: Deleting the interferon receptor (IFNAR1) decreases HSC and MPP progenitor cells, suggesting interferon signaling is required to expand these populations during pregnancy.")

# Run the analysis
analyze_hematopoiesis_data()
