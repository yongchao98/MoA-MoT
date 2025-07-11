import math

# --- Part 1: Calculate the initial proportion of H3K4me3 sites ---

# Given values for Part 1
final_proportion = 11.04  # in percent
decay_rate = 0.10        # 10% per hour
time = 10                # hours

# The formula for exponential decay is: P(t) = P(0) * (1 - r)^t
# We need to find P(0), so we rearrange: P(0) = P(t) / (1 - r)^t
initial_proportion = final_proportion / ((1 - decay_rate) ** time)

print("--- Part 1: Initial Percentage of Trimethylated Sites ---")
print("The relationship between initial and final proportions is described by:")
print(f"Final Proportion = Initial Proportion * (1 - Decay Rate) ^ Time")
print("Plugging in the known values, the equation is:")
print(f"{final_proportion}% = Initial Proportion * (1 - {decay_rate}) ^ {time}")
print("\nTo find the initial proportion, we rearrange the formula:")
print(f"Initial Proportion = {final_proportion}% / (1 - {decay_rate}) ^ {time}")
print(f"Initial Proportion = {final_proportion} / {1 - decay_rate} ^ {time}")
print(f"Initial Proportion = {final_proportion} / {(1 - decay_rate)**time:.6f}")
print(f"\nThe initial percentage of H3K4me3 sites was: {initial_proportion:.2f}%\n")


# --- Part 2: Calculate the impact on gene expression ---

# Given values for Part 2
initial_expression = 200  # RPKM
proportion_decrease = 0.10 # 10% relative decrease

# With a linear relationship, a 10% decrease in the cause (methylation)
# results in a 10% decrease in the effect (expression).
new_expression = initial_expression * (1 - proportion_decrease)

print("--- Part 2: Impact on Gene Expression ---")
print("Assuming a linear relationship, a 10% decrease in methylation proportion leads to a 10% decrease in expression.")
print("The formula to calculate the new expression level is:")
print("New Expression = Initial Expression * (1 - Proportional Decrease)")
print("\nPlugging in the known values, the equation is:")
print(f"New Expression = {initial_expression} RPKM * (1 - {proportion_decrease})")
print(f"New Expression = {initial_expression} RPKM * {1 - proportion_decrease}")
print(f"\nThe new gene expression level is: {new_expression:.0f} RPKM")

# --- Final Answer ---
final_answer = f"The initial percentage of trimethylated sites was {initial_proportion:.2f}%. The new gene expression level is {new_expression:.0f} RPKM."
print(f"\n<<<{final_answer}>>>")