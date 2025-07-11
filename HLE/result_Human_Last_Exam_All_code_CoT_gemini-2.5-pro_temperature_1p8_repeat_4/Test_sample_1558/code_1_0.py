# All concentrations are per ul of blood

# Experiment 1: Effect of Reverse Transcriptase Inhibitors (RTI)
# A positive difference indicates that the treatment caused a decrease compared to the control.
preg_control_rbc_exp1 = 10 * 10**6
preg_rti_rbc_exp1 = 8 * 10**6
rbc_change_rti = preg_control_rbc_exp1 - preg_rti_rbc_exp1

preg_control_bm = 50  # in percent
preg_rti_bm = 30  # in percent
bm_change_rti = preg_control_bm - preg_rti_bm

# Experiment 2: Effect of STING deletion
preg_control_rbc_exp2 = 13 * 10**6
preg_dsting_rbc_exp2 = 8 * 10**6
rbc_change_sting = preg_control_rbc_exp2 - preg_dsting_rbc_exp2

print("--- Analysis of Experiment 1: RTI Effect on Pregnant Mice ---")
print(f"RBC count in control group: {int(preg_control_rbc_exp1)}")
print(f"RBC count in RTI-treated group: {int(preg_rti_rbc_exp1)}")
print("Equation for change in RBC count with RTI treatment:")
print(f"{int(preg_control_rbc_exp1)} (control) - {int(preg_rti_rbc_exp1)} (RTI) = {int(rbc_change_rti)}")
print(f"This decrease of {int(rbc_change_rti)} cells/ul suggests that transposable element activity supports red blood cell production.")
print("-" * 20)

print("--- Analysis of Experiment 2: STING Deletion Effect on Pregnant Mice ---")
print(f"RBC count in control group: {int(preg_control_rbc_exp2)}")
print(f"RBC count in delta STING group: {int(preg_dsting_rbc_exp2)}")
print("Equation for change in RBC count with STING deletion:")
print(f"{int(preg_control_rbc_exp2)} (control) - {int(preg_dsting_rbc_exp2)} (delta STING) = {int(rbc_change_sting)}")
print(f"This decrease of {int(rbc_change_sting)} cells/ul suggests the immune system's STING pathway is essential for enhanced erythropoiesis during pregnancy.")
print("-" * 20)

print("\nConclusion: The data shows that TE activity promotes erythropoiesis via the STING-interferon immune pathway.")
print("The only statement not contradicted by the data is that inducing TEs might be a strategy to combat anemia (low RBCs).")
