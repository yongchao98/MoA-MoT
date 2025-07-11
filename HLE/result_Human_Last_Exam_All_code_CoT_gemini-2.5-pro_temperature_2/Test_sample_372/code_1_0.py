import math

# --- Part 1: Calculate the initial percentage of trimethylated sites ---

# Define the known variables for the exponential decay model
final_percentage_part1 = 11.04  # The percentage of H3K4me3 sites after 10 hours
decay_rate = 0.10             # The rate of turnover per hour (10%)
time = 10                     # The duration in hours

# The formula for exponential decay is: Final = Initial * (1 - rate)^time
# To find the initial percentage, we rearrange the formula: Initial = Final / (1 - rate)^time
initial_percentage = final_percentage_part1 / ((1 - decay_rate)**time)

print("--- Part 1: Calculating the Initial Percentage of Trimethylated Sites ---")
print("This is an exponential decay problem. The formula is: P(t) = P0 * (1 - r)^t")
print("We can solve for the initial percentage (P0) by rearranging the formula: P0 = P(t) / (1 - r)^t")
print(f"\nGiven values:")
print(f"  Final percentage P({time}) = {final_percentage_part1}%")
print(f"  Decay rate r = {decay_rate*100}% per hour")
print(f"  Time t = {time} hours")
print("\nCalculation:")
print(f"  P0 = {final_percentage_part1} / (1 - {decay_rate})^{time}")
print(f"  P0 = {final_percentage_part1} / ({1 - decay_rate})^{time}")
print(f"  P0 = {final_percentage_part1} / {math.pow(1 - decay_rate, time):.4f}")
print(f"Result: The initial percentage of trimethylated sites was {initial_percentage:.2f}%.")


# --- Part 2: Determine the impact on gene expression ---
print("\n" + "---" * 20 + "\n")

# Define the known variables for the gene expression calculation
initial_expression = 200      # Average expression level in RPKM
methylation_decrease = 0.10   # The total decrease in methylation proportion over the period (10%)

# The relationship is linear, so a % decrease in methylation causes the same % decrease in expression.
# New Expression = Initial Expression * (1 - decrease)
final_expression = initial_expression * (1 - methylation_decrease)

print("--- Part 2: Determining the Impact on Gene Expression ---")
print("Assuming a linear relationship between methylation and gene expression, a percentage decrease in H3K4me3 proportion will cause an equal percentage decrease in expression.")
print("\nGiven values:")
print(f"  Initial gene expression = {initial_expression} RPKM")
print(f"  Decrease in H3K4me3 proportion = {methylation_decrease*100}%")
print("\nCalculation:")
print(f"  New Expression = Initial Expression * (1 - Percentage Decrease)")
print(f"  New Expression = {initial_expression} * (1 - {methylation_decrease})")
print(f"  New Expression = {initial_expression} * {1 - methylation_decrease}")
print(f"Result: The new target gene expression level is {final_expression:.0f} RPKM.")

<<<180>>>