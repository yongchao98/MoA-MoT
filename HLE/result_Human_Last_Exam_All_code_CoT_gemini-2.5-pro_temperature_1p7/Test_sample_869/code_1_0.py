import math

# Step 1: Define the parameters
N = 25  # Total number of items
T = 5   # Number of types of items, also the number of individuals
C = 5   # Number of copies of each type of item

# Step 2: Calculate the total number of possible distributions (S)
# S is the number of permutations of a multiset of 25 items,
# where there are 5 identical items of 5 different types.
# S = 25! / (5! * 5! * 5! * 5! * 5!) = 25! / (5!)^5
try:
    s_denominator = math.factorial(C)**T
    S = math.factorial(N) // s_denominator
except OverflowError:
    # Use floating point for very large numbers if precise integer is not needed
    # but for this problem, we need precision.
    # Python's integers handle arbitrary size, so this is unlikely unless N is huge.
    # For calculation, it's better to use logs for intermediate steps, but for the final
    # values of F and S themselves, we can compute them directly.
    s_denominator = float(math.factorial(C))**T
    S = float(math.factorial(N)) / s_denominator


# Step 3: Calculate the number of favorable distributions (F)
# Based on the assumption that a favorable outcome occurs only when each individual
# receives all 5 items of one unique type.
# The number of ways to assign the 5 unique types to the 5 individuals is 5!.
# For each assignment, there's only 1 way the items can be arranged.
# So, F = 5!
F = math.factorial(T)

# Step 4: Calculate the probability P = F / S
P = F / S

# Step 5: Output the results
# The problem asks to output the final equation with each number.
print("Calculating the probability P = F / S")
print(f"Total number of distributions (S) = {N}! / ({C}!)**{T} = {S}")
print(f"Number of favorable distributions (F) = {T}! = {F}")
print(f"The final equation is: P = {F} / {S}")
print(f"The calculated probability P is: {P}")
