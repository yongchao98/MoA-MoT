# Data from Experiment 1 for Red Blood Cells in pregnant mice.
# All values are the coefficients for the factor 10^6 per ul.
pregnant_control_rbc = 10
pregnant_rti_rbc = 8

# Calculate the percentage decrease in Red Blood Cells.
# The formula is ((initial - final) / initial) * 100
decrease = pregnant_control_rbc - pregnant_rti_rbc
percentage_decrease = (decrease / pregnant_control_rbc) * 100

# Print the result in a sentence that includes all the numbers from the equation.
print("From the data in experiment 1, we can calculate the percentage decrease in Red Blood Cells in pregnant mice treated with RTI.")
print(f"The equation is: (({pregnant_control_rbc} - {pregnant_rti_rbc}) / {pregnant_control_rbc}) * 100")
print(f"The resulting decrease is {percentage_decrease:.1f}%.")
