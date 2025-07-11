import math

# Part 1: Calculate the initial percentage of trimethylated sites

# Known values
final_proportion_percent = 11.04
time_hours = 10
decay_rate_per_hour_percent = 10

# Convert percentages to decimals for calculation
final_proportion = final_proportion_percent / 100
decay_rate = decay_rate_per_hour_percent / 100

# The exponential decay formula is P(t) = P(0) * e^(-k*t)
# We rearrange it to solve for the initial proportion P(0): P(0) = P(t) / e^(-k*t)
initial_proportion = final_proportion / math.exp(-decay_rate * time_hours)

# Convert the result back to a percentage
initial_percentage = initial_proportion * 100

print("Part 1: Initial Trimethylation Percentage")
print(f"The equation for the initial proportion is: P(0) = {final_proportion_percent}% / e^(-{decay_rate_per_hour_percent}% * {time_hours} hours)")
print(f"The calculated percentage of sites trimethylated at the start is: {initial_percentage:.2f}%\n")


# Part 2: Determine the impact on gene expression

# Known values
initial_expression_rpkm = 200
proportion_decrease_percent = 10

# Calculate the remaining proportion factor
proportion_multiplier = (100 - proportion_decrease_percent) / 100

# Assuming a linear relationship, the new expression level is the initial level times the new proportion factor
new_expression = initial_expression_rpkm * proportion_multiplier

print("Part 2: Impact on Gene Expression")
print(f"The equation for the new expression level is: New Expression = {initial_expression_rpkm} RPKM * (1 - {proportion_decrease_percent / 100})")
print(f"With a {proportion_decrease_percent}% decrease in H3K4me3 sites, the new gene expression level is: {new_expression:.2f} RPKM")

# The final answer is the initial percentage from Part 1.
# The user might be looking for a single value. Part 1 asks "what percentage of sites is trimethylated at the start", which is the primary question.
# The second part is a subsequent calculation. I'll use the result from the first part as the answer.
final_answer = initial_percentage
print(f"<<<{final_answer:.2f}>>>")