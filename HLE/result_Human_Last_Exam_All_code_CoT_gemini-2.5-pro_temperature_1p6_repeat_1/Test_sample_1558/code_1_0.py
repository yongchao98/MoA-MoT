# The plan is to calculate the percentage change in red blood cells (RBCs)
# in pregnant mice when treated with Reverse Transcriptase Inhibitors (RTI),
# based on the data from the first experiment. This calculation will highlight
# the impact of inhibiting transposable elements on erythropoiesis (RBC production)
# and help validate the chosen conclusion.

# Data from Experiment 1 for pregnant mice:
# Red Blood Cell count in the control group.
pregnant_rbc_control = 10.0  # value is 10 x 10^6 per ul

# Red Blood Cell count in the RTI-treated group.
pregnant_rbc_rti = 8.0  # value is 8 x 10^6 per ul

# The formula for percentage change is: ((new_value - old_value) / old_value) * 100
percentage_change = ((pregnant_rbc_rti - pregnant_rbc_control) / pregnant_rbc_control) * 100

print("This script calculates the effect of inhibiting transposable elements on red blood cell (RBC) counts in pregnant mice.")
print(f"The control group of pregnant mice had an RBC count of {int(pregnant_rbc_control)}x10^6 per ul.")
print(f"The RTI-treated group of pregnant mice had an RBC count of {int(pregnant_rbc_rti)}x10^6 per ul.")
print("\nTo find the percentage change, we use the equation: ((Treated Value - Control Value) / Control Value) * 100%")
print("\nFinal Equation using the experimental data:")
# The prompt requires outputting each number in the final equation.
print(f"(({pregnant_rbc_rti} - {pregnant_rbc_control}) / {pregnant_rbc_control}) * 100 = {percentage_change:.1f}%")

print("\nThe calculation shows a 20.0% decrease in RBCs in pregnant mice when reverse transcriptase (and thus transposable element activity) is inhibited.")
print("This suggests that the natural activity of transposable elements during pregnancy helps support red blood cell production, counteracting anemia. Therefore, the induction of transposons might be a potential therapeutic approach for treating anemia.")
