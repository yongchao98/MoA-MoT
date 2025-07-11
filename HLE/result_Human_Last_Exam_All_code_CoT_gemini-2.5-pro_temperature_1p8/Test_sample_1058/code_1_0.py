import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("Analyzing the experimental data to find the correct answer:\n")

# Step 1: Analyze neuronal activity from Population Spike (PS) amplitudes
ps_alcohol_preferring = -0.38  # mV
ps_sucrose_preferring = -0.17  # mV

print("Step 1: Comparing Neuronal Activity")
print(f"PS amplitude in alcohol-preferring rats: {ps_alcohol_preferring} mV")
print(f"PS amplitude in sucrose-preferring rats: {ps_sucrose_preferring} mV")
if abs(ps_alcohol_preferring) > abs(ps_sucrose_preferring):
    activity_conclusion = "increased neuronal activity/excitability"
else:
    activity_conclusion = "decreased neuronal activity/excitability"
print(f"Conclusion: A larger PS amplitude magnitude indicates {activity_conclusion}.\n")


# Step 2: Analyze the Slc6a11 knockdown experiment
ps_shRNA_knockdown = -0.37 # mV
ps_control = -0.16 # mV
print("Step 2: Analyzing the Gene Knockdown Experiment")
print("The gene Slc6a11 codes for a GABA transporter (GAT), which removes GABA from the extracellular space.")
print("Knocking down this gene using shRNA results in a PS amplitude of {} mV, mimicking the alcohol-preferring rats ({} mV).".format(ps_shRNA_knockdown, ps_alcohol_preferring))
gaba_conclusion = "higher in the alcohol-preferring rats, as reduced transporter function leads to less GABA removal."
print(f"Conclusion: This implies the level of extracellular GABA is {gaba_conclusion}\n")


# Step 3: Analyze tonic inhibition
print("Step 3: Deducing Tonic Inhibition")
print("Tonic inhibition is driven by extracellular GABA levels.")
tonic_inhibition_conclusion = "increased tonic inhibition in alcohol-preferring rats."
print(f"Conclusion: Higher extracellular GABA leads to {tonic_inhibition_conclusion}\n")


# Step 4: Proposing a therapeutic strategy
print("Step 4: Proposing a Therapeutic Strategy")
print("The issue appears to be excessive GABA signaling leading to network hyperexcitability (disinhibition).")
therapy_conclusion = "A GABA receptor ANTAGONIST would be needed to block the excess GABA effect. An AGONIST would likely make it worse."
print(f"Conclusion: {therapy_conclusion}\n")


# Step 5: Evaluating option G
print("Step 5: Evaluating the Answer Choices")
print("Option G states:")
print("  - Alcohol-preferring rats show increased neuronal activity in the amygdala. (Correct)")
print("  - GABA receptor antagonists may decrease the self-administration of alcohol. (Correct)")
print("  - The level of extracellular GABA is higher in the alcohol-preferring rats. (Correct)")
print("  - Alcohol-preferring rats have increased tonic inhibition in the amygdala. (Correct)")
print("\nAll statements in option G are supported by the analysis of the experimental data.")


# Restore stdout
sys.stdout = old_stdout
# Print the captured output
print(captured_output.getvalue())
print("<<<G>>>")