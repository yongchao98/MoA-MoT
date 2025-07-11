import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Given values
p_final_percent = 11.04  # Final percentage after 10 hours
decay_rate_per_hour = 0.10 # 10% per hour
time_hours = 10          # 10 hours

# The formula for exponential decay is P(t) = P0 * (1 - r)^t
# We rearrange it to solve for the initial percentage, P0 = P(t) / (1 - r)^t

# Calculate the decay factor
decay_factor = (1 - decay_rate_per_hour) ** time_hours

# Calculate the initial percentage
initial_percentage = p_final_percent / decay_factor

print("--- Part 1: Initial Methylation Percentage ---")
print("The calculation for the initial percentage (P0) is based on the final percentage (P(t)), the decay rate (r), and time (t).")
print(f"P0 = P(t) / (1 - r)^t")
print(f"P0 = {p_final_percent} / (1 - {decay_rate_per_hour})^{time_hours}")
print(f"The initial percentage of H3K4me3 sites was: {initial_percentage:.2f}%\n")


# Part 2: Determine the impact on gene expression

# Given values
initial_expression = 200  # RPKM
proportion_decrease = 0.10 # 10% decrease

# Assuming a linear relationship, a 10% decrease in the methylation proportion
# leads to a 10% decrease in the gene expression level.
final_expression = initial_expression * (1 - proportion_decrease)

print("--- Part 2: Impact on Gene Expression ---")
print("The calculation for the new gene expression level is based on the initial level and the percentage decrease.")
print(f"New Expression = Initial Expression * (1 - Decrease)")
print(f"New Expression = {initial_expression} * (1 - {proportion_decrease})")
print(f"The new gene expression level would be: {final_expression:.2f} RPKM")

# Final combined answer
final_answer_part1 = round(initial_percentage, 2)
final_answer_part2 = round(final_expression, 2)
# The final answer format is not directly printed to the console in this code block
# but is provided after the code block as requested.
