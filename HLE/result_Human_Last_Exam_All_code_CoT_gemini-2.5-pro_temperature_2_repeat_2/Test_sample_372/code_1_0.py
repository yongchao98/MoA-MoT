import math

# Part 1: Initial methylation percentage
print("Part 1: Initial Trimethylation Percentage")
print("-" * 40)

# Given values for Part 1
rate_turnover = 0.10  # 10% per hour
time_hours = 10
final_proportion_me3 = 0.1104 # 11.04%

# The decay model is P_final = P_initial * (1 - rate)^time
# We solve for P_initial = P_final / (1 - rate)^time
initial_proportion_me3 = final_proportion_me3 / ((1 - rate_turnover)**time_hours)
initial_percentage_me3 = initial_proportion_me3 * 100

print("The formula to find the initial proportion (P_initial) is:")
print("P_initial = P_final / (1 - rate)^time\n")
print("Calculation:")
print(f"P_initial = {final_proportion_me3} / (1 - {rate_turnover})^{time_hours}")
print(f"P_initial = {final_proportion_me3} / ({(1 - rate_turnover)**time_hours})")
print(f"P_initial = {initial_proportion_me3:.4f}")
print(f"\nThe initial percentage of H3K4me3 sites was {initial_percentage_me3:.2f}%.")
print("\n" + "=" * 40 + "\n")


# Part 2: Impact on gene expression
print("Part 2: Impact on Gene Expression")
print("-" * 40)

# Given values for Part 2
initial_expression = 200 # RPKM
decrease_factor = 0.10 # 10% relative decrease

# With a linear relationship (Expression = k * Proportion), a 10% relative decrease
# in proportion results in a 10% relative decrease in expression.
# new_expression = initial_expression * (1 - 0.10)
new_expression = initial_expression * (1 - decrease_factor)

print("Assuming a linear relationship, the new expression level (E_new) is calculated as:")
print("E_new = E_initial * (1 - decrease_percentage)\n")
print("Calculation:")
print(f"E_new = {initial_expression} RPKM * (1 - {decrease_factor})")
print(f"E_new = {initial_expression} RPKM * {1 - decrease_factor}")
print(f"E_new = {new_expression:.2f} RPKM")
print(f"\nThe new average gene expression level would be {new_expression:.2f} RPKM.")

<<<180.00>>>