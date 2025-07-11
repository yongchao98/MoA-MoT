import math

# Part 1: Calculate the initial percentage of H3K4me3 sites.
# The turnover of H3K4me3 is modeled by the exponential decay formula: P(t) = P(0) * e^(-k*t)
# We are given P(10) = 11.04%, k = 10%/hr = 0.10, and t = 10 hr.
# We can find the initial proportion P(0) by rearranging the formula: P(0) = P(t) * e^(k*t).

print("### Part 1: Calculating the Initial H3K4me3 Percentage ###")

# Given values for Part 1
final_proportion_pct = 11.04
decay_rate = 0.10  # 10% per hour as a decimal
time = 10          # hours

# Convert final proportion from percentage to a decimal for the calculation
final_proportion = final_proportion_pct / 100.0

# Calculate the initial proportion using the rearranged formula
initial_proportion = final_proportion * math.exp(decay_rate * time)

# Convert the initial proportion back to a percentage for the output
initial_proportion_pct = initial_proportion * 100.0

print(f"The calculation for the initial proportion P(0) is:")
print(f"P(0) = P({time}) * e^(k * t)")
print(f"P(0) = {final_proportion_pct}% * e^({decay_rate} * {time})")
print(f"The initial percentage of trimethylated sites is {initial_proportion_pct:.2f}%.")
print("-" * 30)

# Part 2: Determine the impact on gene expression.
# The relationship between proportion and expression is linear. We'll find the new expression
# level after a 10 percentage point decrease in H3K4me3.

print("\n### Part 2: Determining the Impact on Gene Expression ###")

# Given values for Part 2
initial_expression = 200      # RPKM
proportion_decrease_pct = 10  # A decrease of 10 percentage points

# Calculate the new proportion after the decrease
new_proportion_pct = initial_proportion_pct - proportion_decrease_pct
new_proportion = new_proportion_pct / 100.0

# Assume a linear relationship: Expression = m * Proportion.
# First, calculate the proportionality constant 'm' from the initial state.
m = initial_expression / initial_proportion

# Then, calculate the new expression level using the new proportion.
new_expression = m * new_proportion

# Finally, calculate the impact, which is the change in expression.
expression_change = new_expression - initial_expression

print(f"The initial proportion of H3K4me3 is {initial_proportion_pct:.2f}% with an expression of {initial_expression} RPKM.")
print(f"A decrease by {proportion_decrease_pct} percentage points results in a new proportion of {new_proportion_pct:.2f}%.")
print(f"\nAssuming a linear model (Expression = m * Proportion), the new expression level is calculated as:")
print(f"New Expression = (Initial Expression / Initial Proportion) * New Proportion")
print(f"New Expression = ({initial_expression} RPKM / {initial_proportion_pct:.2f}%) * {new_proportion_pct:.2f}%")
print(f"New Expression = {new_expression:.2f} RPKM")
print(f"\nThe impact is the change in expression:")
print(f"Impact = New Expression - Initial Expression")
print(f"Impact = {new_expression:.2f} RPKM - {initial_expression} RPKM")
print(f"Impact = {expression_change:.2f} RPKM")

<<< -66.67 >>>