import math

# Part 1: Calculate the initial percentage of H3K4me3 sites

# Given values
rate_turnover = 0.10  # 10% per hour
time_hours = 10.0      # 10 hours
final_proportion = 0.1104  # 11.04%

print("--- Part 1: Calculating the Initial Percentage of Trimethylated Sites ---")
print("The model for histone turnover follows the exponential decay formula: P(t) = P(0) * e^(-r*t)")
print("We need to solve for the initial proportion, P(0).")
print(f"P({int(time_hours)}) = P(0) * e^(-{rate_turnover} * {int(time_hours)})")
print(f"{final_proportion} = P(0) * e^(-{rate_turnover * time_hours})")

# To find P(0), we rearrange the formula: P(0) = P(t) * e^(r*t)
initial_proportion = final_proportion * math.exp(rate_turnover * time_hours)
initial_percentage = initial_proportion * 100

print(f"P(0) = {final_proportion} * e^({rate_turnover} * {int(time_hours)})")
print(f"P(0) = {initial_proportion:.4f}")
print(f"\nThe percentage of sites trimethylated at the start was: {initial_percentage:.2f}%\n")


# Part 2: Determine the impact on gene expression

# Given values
initial_expression = 200  # RPKM
proportion_decrease = 10  # 10%

print("--- Part 2: Determining the Impact on Gene Expression ---")
print("Assuming a linear relationship between methylation proportion and gene expression.")
print("A 10% decrease in H3K4me3 proportion will cause a 10% decrease in expression.")

# Calculate the new expression level
final_expression = initial_expression * (1 - proportion_decrease / 100)

print(f"New Expression = Initial Expression * (1 - Percentage Decrease / 100)")
print(f"New Expression = {initial_expression} * (1 - {proportion_decrease} / 100)")
print(f"New Expression = {initial_expression} * {1 - proportion_decrease / 100}")
print(f"New Expression = {final_expression:.2f} RPKM")

print("\nThe impact on gene expression is a decrease from 200 RPKM to 180 RPKM.")
