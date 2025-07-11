import math

# Part 1: Calculate the initial percentage of H3K4me3 sites.

# Given values
p_t = 0.1104  # Proportion of H3K4me3 after 10 hours
rate = 0.10   # Rate of turnover per hour
time = 10     # Time in hours

# We use the formula for exponential decay: P(t) = P(0) * (1 - r)^t
# We need to find P(0), the initial proportion.
# Rearranging the formula: P(0) = P(t) / (1 - r)^t

one_minus_rate_to_t = math.pow(1 - rate, time)
p_0 = p_t / one_minus_rate_to_t

# Part 2: Calculate the impact on gene expression.

# Given values
initial_expression = 200  # RPKM
proportion_decrease = 0.10 # 10% decrease

# Assuming a linear relationship, a 10% decrease in the methylation proportion
# leads to a 10% decrease in gene expression.
final_expression = initial_expression * (1 - proportion_decrease)

# --- Output ---
print("Part 1: Calculating the initial percentage of trimethylated sites")
print("Formula: Initial_Percentage = Final_Percentage / (1 - Rate)^Time")
print(f"Calculation: Initial_Percentage = {p_t:.4f} / (1 - {rate:.2f})^{time}")
print(f"Resulting initial proportion: {p_0:.4f}, which is {p_0 * 100:.2f}%")
print("-" * 50)
print("Part 2: Determining the impact on gene expression")
print("A 10% decrease in methylation proportion causes a 10% decrease in expression.")
print("Formula: New_Expression = Initial_Expression * (1 - Decrease_Percentage)")
print(f"Calculation: New_Expression = {initial_expression} * (1 - {proportion_decrease:.2f})")
print(f"Resulting new expression level: {final_expression:.2f} RPKM")

# Format the final answer as requested
final_answer_str = f"The initial percentage of trimethylated sites is {p_0 * 100:.2f}%. The target gene expression would decrease to {final_expression:.2f} RPKM."
print(f"\n<<<The initial percentage of trimethylated sites is {p_0 * 100:.2f}%. The target gene expression would decrease to {final_expression:.2f} RPKM.>>>")