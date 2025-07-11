import math

# --- Data from the experiments ---

# Experiment 1: RTI treatment
exp1_rbc = {
    "non_pregnant_control": 13 * 10**6,
    "non_pregnant_rti": 13 * 10**6,
    "pregnant_control": 10 * 10**6,
    "pregnant_rti": 8 * 10**6,
}

# Experiment 3: ifnar1 deletion (progenitor cells in spleen)
exp3_hsc = {
    "pregnant_control": 0.003, # percent
    "pregnant_delta_ifnar1": 0.002, # percent
}


# --- Analysis ---

print("Step 1: Analyzing the effect of TE inhibition on Red Blood Cells in pregnant mice (Experiment 1).")
preg_control_rbc = exp1_rbc["pregnant_control"]
preg_rti_rbc = exp1_rbc["pregnant_rti"]
print(f"The number of Red Blood Cells in pregnant mice decreased from {preg_control_rbc:.0e} per ul (control) to {preg_rti_rbc:.0e} per ul (RTI treated).")
print("Conclusion 1: Since inhibiting transposable elements (TEs) with RTI caused a drop in Red Blood Cells, this suggests that the normal activity of TEs supports or increases the number of Red Blood Cells in pregnant mice.")
print("-" * 20)


print("Step 2: Analyzing the effect of Interferon signaling on blood cell precursors (Experiment 3).")
preg_control_hsc = exp3_hsc["pregnant_control"]
preg_ifnar1_hsc = exp3_hsc["pregnant_delta_ifnar1"]
print(f"The percentage of HSC precursor cells in the spleen of pregnant mice decreased from {preg_control_hsc}% (control) to {preg_ifnar1_hsc}% (delta ifnar1).")
print("Conclusion 2: Blocking the Interferon pathway reduces the number of blood cell precursors. This suggests the Interferon pathway is important for hematopoiesis.")
print("-" * 20)


print("Step 3: Evaluating Answer Choice A.")
print("Answer A states: 'Increased activity of transposable elements increases the number of red blood cells in pregnant mice. Interferon does not increase the number of red blood cells in pregnant mice.'")
print("\nEvaluation of the first part: 'Increased activity of transposable elements increases the number of red blood cells in pregnant mice.'")
print("Our analysis of Experiment 1 supports this. Inhibiting TEs with RTI lowered the red blood cell count from 1.00e+07 to 8.00e+06.")

print("\nEvaluation of the second part: 'Interferon does not increase the number of red blood cells in pregnant mice.'")
print("Our analysis of Experiment 3 shows that Interferon signaling is needed for precursor cell production. While a lack of precursors should logically lead to fewer red blood cells, Experiment 3 did not directly measure the final number of red blood cells. Based *only* on the data provided, there is no experiment showing Interferon increasing red blood cells. Therefore, this statement is not contradicted by the direct evidence in the text.")

print("\nFinal verdict: Answer A is the most accurate conclusion based on the provided experimental results.")
