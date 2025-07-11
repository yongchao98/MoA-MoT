import math

# Part 1: Calculate the initial percentage of trimethylated sites.
print("--- Part 1: Initial Percentage of H3K4me3 Sites ---")

# Given values
final_proportion_p1 = 11.04  # in percent
rate = 0.10  # 10% per hour
time = 10  # in hours

# The process is modeled using the continuous exponential decay formula: A(t) = A(0) * e^(-k*t)
# We need to solve for the initial amount, A(0).
# A(0) = A(t) / e^(-k*t) = A(t) * e^(k*t)

initial_proportion_p1 = final_proportion_p1 * math.exp(rate * time)

print(f"The final proportion of H3K4me3 after {time} hours is {final_proportion_p1}%.")
print(f"The decay rate is {rate*100}% per hour.")
print("Using the formula A(0) = A(t) * e^(k*t):")
print(f"Initial Percentage = {final_proportion_p1} * e^({rate} * {time})")
print(f"Initial Percentage = {initial_proportion_p1:.2f}%")
print(f"Therefore, the initial percentage of trimethylated sites was approximately {round(initial_proportion_p1)}%.")

print("\n" + "="*50 + "\n")

# Part 2: Determine the impact on gene expression.
print("--- Part 2: Impact on Target Gene Expression ---")

# Given values
initial_expression = 200  # in RPKM
proportion_decrease = 0.10  # 10% decrease

# Due to the linear relationship, a 10% decrease in methylation proportion
# leads to a 10% decrease in gene expression.
final_expression = initial_expression * (1 - proportion_decrease)

print(f"The initial average gene expression is {initial_expression} RPKM.")
print(f"The proportion of H3K4me3 sites decreases by {proportion_decrease*100}%.")
print("Assuming a linear relationship, the new expression level is calculated as:")
print(f"New Expression = Initial Expression * (1 - Proportion Decrease)")
print(f"New Expression = {initial_expression} * (1 - {proportion_decrease}) = {int(final_expression)} RPKM")
print(f"The impact is a reduction of the target gene expression to {int(final_expression)} RPKM.")

<<<180>>>