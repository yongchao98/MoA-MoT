import math

# Part 1: Calculate the initial percentage of H3K4me3

# Given values for Part 1
P_final_percent = 11.04  # Percentage of H3K4me3 after 10 hours
rate = 0.10              # 10% turnover rate per hour
time = 10                # 10 hours

# The formula for continuous decay is P(t) = P0 * e^(-r*t)
# We rearrange it to find the initial percentage P0 = P(t) / e^(-r*t)
P_initial = P_final_percent / math.exp(-rate * time)

print("--- Part 1: Initial Trimethylation Percentage ---")
print(f"The equation to find the initial percentage (P0) is: P0 = {P_final_percent}% / e^(-{rate} * {time})")
print(f"Calculation: P0 = {P_final_percent} / {math.exp(-rate * time):.4f}")
print(f"The initial percentage of trimethylated sites (P0) was: {P_initial:.2f}%\n")


# Part 2: Determine the impact on gene expression

# Given values for Part 2
initial_expression = 200  # RPKM
methylation_decrease_percent = 10 # 10% decrease in the proportion of sites

# Assuming a linear relationship, a 10% decrease in methylation proportion
# leads to a 10% decrease in gene expression.
new_expression = initial_expression * (1 - methylation_decrease_percent / 100)

print("--- Part 2: Impact on Gene Expression ---")
print(f"The equation to find the new expression level (E_new) is: E_new = {initial_expression} * (1 - {methylation_decrease_percent} / 100)")
print(f"Calculation: E_new = {initial_expression} * {1 - methylation_decrease_percent / 100}")
print(f"With a {methylation_decrease_percent}% decrease in H3K4me3 sites, the new gene expression level is: {new_expression:.2f} RPKM")

# The final answer required is the result of the second part of the question.
final_answer = new_expression
# print(f"<<<{final_answer}>>>")