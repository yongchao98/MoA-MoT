#
# Plan:
# 1. Define variables for the Red Blood Cell (RBC) counts in pregnant mice from Experiment 1.
#    - One variable for the control group.
#    - One variable for the group treated with Reverse Transcriptase Inhibitors (RTI).
# 2. Calculate the difference in RBC counts between the control and treated groups. This difference represents the approximate number of RBCs supported by transposable element activity.
# 3. Print the calculation in a clear, human-readable format, showing all the numbers involved in the equation. This calculation supports the conclusion that transposon activity is beneficial for RBC production in this context.
#

# Step 1: Define RBC counts for pregnant mice (per microliter) from Experiment 1
rbc_pregnant_control = 10 * 10**6
rbc_pregnant_rti = 8 * 10**6

# Step 2: Calculate the reduction in RBCs due to RTI treatment
rbc_reduction = rbc_pregnant_control - rbc_pregnant_rti

# Step 3: Print the full equation and the result
print("This script calculates the reduction in Red Blood Cells (RBCs) in pregnant mice when treated with RTI, based on the data from Experiment 1.")
print("This shows that the activity blocked by RTI was responsible for maintaining a higher RBC count.")
print("\n--- Calculation ---")
print(f"RBCs in Pregnant Control Mice: {rbc_pregnant_control:,}")
print(f"RBCs in Pregnant Mice with RTI: {rbc_pregnant_rti:,}")
print(f"Equation for the reduction: {rbc_pregnant_control:,} - {rbc_pregnant_rti:,} = {rbc_reduction:,}")
print(f"\nResult: The presence of transposable element activity accounts for an additional {rbc_reduction:,} red blood cells per microliter in pregnant mice, helping to counteract anemia.")
