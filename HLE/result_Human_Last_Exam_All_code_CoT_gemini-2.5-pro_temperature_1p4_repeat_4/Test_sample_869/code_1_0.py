import math

# Step 1: Calculate the total number of possible distributions (S)
# S = 25! / (5! * 5! * 5! * 5! * 5!)
s_numerator = math.factorial(25)
s_denominator = math.factorial(5)**5
S = s_numerator / s_denominator

# Step 2: Calculate the number of favorable distributions (F)
# A favorable distribution is one where each person i gets all 5 items of a unique type pi(i).
# There are 5! ways to assign a unique type to each person.
# Each specific assignment corresponds to exactly one sequence of items.
# So, F = 5!
F = math.factorial(5)

# Step 3: Calculate the probability P = F / S
P = F / S

# Step 4: Print the results in the specified format.
# We cast the large numbers to integers for clean printing.
print(f"Total number of distributions (S): {int(S)}")
print(f"Number of favorable distributions (F): {int(F)}")
print(f"Probability (P) = F / S")
print(f"P = {int(F)} / {int(S)}")
print(f"P = {P}")