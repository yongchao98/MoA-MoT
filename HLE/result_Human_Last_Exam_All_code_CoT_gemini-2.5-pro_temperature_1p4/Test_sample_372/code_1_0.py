import math

# --- Part 1: Calculate the initial percentage of H3K4me3 ---

print("--- Part 1: Initial Percentage of H3K4me3 ---")

# Known variables
p_final = 11.04  # Final percentage after 10 hours
rate = 0.10      # 10% turnover rate per hour
time = 10        # 10 hours

# The process is modeled by the exponential decay formula: P(t) = P0 * e^(-k*t)
# We need to solve for the initial percentage, P0.
# P0 = P(t) / e^(-k*t)

# Calculate the initial percentage
p_initial = p_final / math.exp(-rate * time)

# Print the explanation and the result
print("The calculation is based on the exponential decay model: P(t) = P0 * e^(-k*t)")
print(f"Given P({time}) = {p_final}%, the rate k = {rate}, and time t = {time} hours.")
print(f"The equation to solve for the initial percentage (P0) is:")
print(f"{p_final} = P0 * e^(-{rate} * {time})")
print(f"P0 = {p_final} / e^(-{rate * time})")
print(f"The initial percentage of trimethylated sites (P0) was: {p_initial:.2f}%\n")


# --- Part 2: Calculate the impact on gene expression ---

print("--- Part 2: Impact on Gene Expression ---")

# Known variables
initial_expression = 200  # RPKM
methylation_decrease = 10 # 10% decrease

# Due to the linear relationship, a 10% decrease in methylation causes a 10% decrease in expression.
decrease_factor = methylation_decrease / 100.0

# Calculate the final expression and the change in expression
final_expression = initial_expression * (1 - decrease_factor)
expression_change = initial_expression - final_expression

# Print the explanation and the result
print("Assuming a linear relationship between methylation and gene expression.")
print(f"A decrease in H3K4me3 proportion by {methylation_decrease}% will cause a corresponding decrease in gene expression.")
print(f"Initial gene expression: {initial_expression} RPKM")
print(f"The final expression level is calculated as: {initial_expression} * (1 - {decrease_factor})")
print(f"Final gene expression: {final_expression:.2f} RPKM")
print(f"This represents a decrease of {expression_change:.2f} RPKM in gene expression.")

<<<30.01>>>