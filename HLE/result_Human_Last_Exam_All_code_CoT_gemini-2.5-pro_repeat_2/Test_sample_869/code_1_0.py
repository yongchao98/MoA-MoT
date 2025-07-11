import math

# Step 1: Calculate the total number of distributions (S)
# S = 25! / (5! * 5! * 5! * 5! * 5!)
fact_5 = math.factorial(5)
fact_25 = math.factorial(25)
S = fact_25 / (fact_5 ** 5)

# Step 2: Calculate the number of favorable distributions (F)
# F is the number of ways to assign 5 pure hands of unique types to 5 individuals.
F = fact_5

# Step 3: Calculate the probability P = F / S
P = F / S

# Print the components of the equation as requested.
# Using integer representation for F and S for clarity in the equation.
print(f"F = {int(F)}")
print(f"S = {int(S)}")
print(f"P = F / S = {int(F)} / {int(S)} = {P}")
