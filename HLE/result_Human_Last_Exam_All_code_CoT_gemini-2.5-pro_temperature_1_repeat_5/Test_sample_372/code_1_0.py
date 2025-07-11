import math

# --- Part 1: Calculating the Initial Percentage of H3K4me3 ---

# Known values from the problem
final_proportion_percent = 11.04
decay_rate_per_hour = 0.10
time_hours = 10

# Convert the final percentage to a decimal for calculation
final_proportion_decimal = final_proportion_percent / 100.0

# The formula for exponential decay is A(t) = A(0) * (1 - r)^t.
# We rearrange it to find the initial amount A(0): A(0) = A(t) / (1 - r)^t
initial_proportion_decimal = final_proportion_decimal / ((1 - decay_rate_per_hour) ** time_hours)

# Convert the result back to a percentage
initial_proportion_percent = initial_proportion_decimal * 100

print("--- Part 1: Initial Percentage of Trimethylated Sites ---")
print("To find the initial percentage, we use the formula for exponential decay, A(0) = A(t) / (1 - r)^t.")
print("\nCalculation for the initial percentage:")
# The final code outputs each number in the final equation as requested
print(f"Initial Percentage = {final_proportion_percent}% / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"Initial Percentage = {round(initial_proportion_percent, 2)}%\n")


# --- Part 2: Calculating the Impact on Gene Expression ---

# Known values from the problem
initial_expression_rpkm = 200
proportion_decrease = 0.10 # This represents a 10% decrease

# Due to the linear relationship, a 10% decrease in methylation proportion
# results in a 10% decrease in gene expression.
# New Expression = Initial Expression * (1 - decrease)
new_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

print("--- Part 2: Impact on Gene Expression ---")
print("Given a linear relationship, the new expression level is calculated based on the proportional decrease in methylation.")
print("\nCalculation for the new gene expression level:")
# The final code outputs each number in the final equation as requested
print(f"New Expression (RPKM) = {initial_expression_rpkm} * (1 - {proportion_decrease})")
print(f"New Expression (RPKM) = {round(new_expression_rpkm, 2)}\n")

print(f"Summary:")
print(f"The initial percentage of sites that were trimethylated is {round(initial_proportion_percent, 2)}%.")
print(f"If the proportion of H3K4me3 sites decreases by 10%, the new gene expression level is {round(new_expression_rpkm, 2)} RPKM.")
