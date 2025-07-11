import math

# Plan:
# 1. Define the growth rates phi and alpha from the characteristic equations of the recurrences.
# 2. Calculate the limit L = 2/ln(phi) + 2/ln(alpha).
# 3. Compute the final value 10^4 * L and find its integer part.
# 4. Print the numbers used in the final equation as requested.

# Step 1: Define the growth rates
# For k=5, the characteristic equation is r^2 - 5*r + 1 = 0
phi = (5 + math.sqrt(21)) / 2
# For k=13, the characteristic equation is r^2 - 13*r + 1 = 0
alpha = (13 + math.sqrt(165)) / 2

# Step 2: Calculate the limit L
c1 = 2 / math.log(phi)
c2 = 2 / math.log(alpha)
limit_L = c1 + c2

# Step 3: Compute the final value and its integer part
final_value = 10000 * limit_L
integer_part = int(final_value)

# Step 4: Output the results showing the "final equation"
print("The limit is L = lim F(N)/ln(N) = C1 + C2.")
print(f"The first component of the limit is C1 = 2 / ln((5 + sqrt(21))/2)")
print(f"C1 = {c1}")
print(f"The second component of the limit is C2 = 2 / ln((13 + sqrt(165))/2)")
print(f"C2 = {c2}")
print(f"The total limit is L = {c1} + {c2} = {limit_L}")
print("\nThe final result is the integer part of the expression: 10^4 * L")
print(f"10000 * {limit_L} = {final_value}")
print(f"The integer part is: {integer_part}")
