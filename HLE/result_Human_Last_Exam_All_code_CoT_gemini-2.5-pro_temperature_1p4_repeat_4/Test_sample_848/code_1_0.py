import math

# The problem is to find the integer part of 10^4 * L, where L is the limit.
# The limit L is the sum of contributions from two families of solutions, corresponding to k=13 and k=5.
# The formula for the limit contribution from a family with characteristic root r is 2 / ln(r).

print("--- Calculation for the k=13 family ---")
# The characteristic root r_13 is the larger root of r^2 - 13r + 1 = 0.
r_13 = (13 + 3 * math.sqrt(17)) / 2
print(f"Characteristic root r_13 = (13 + 3*sqrt(17)) / 2 = {r_13}")
# The natural logarithm of r_13.
ln_r_13 = math.log(r_13)
print(f"ln(r_13) = {ln_r_13}")
# The contribution to the limit from this family.
lim_13 = 2 / ln_r_13
print(f"Limit contribution = 2 / {ln_r_13} = {lim_13}")

print("\n--- Calculation for the k=5 family ---")
# The characteristic root r_5 is the larger root of r^2 - 5r + 1 = 0.
r_5 = (5 + math.sqrt(21)) / 2
print(f"Characteristic root r_5 = (5 + sqrt(21)) / 2 = {r_5}")
# The natural logarithm of r_5.
ln_r_5 = math.log(r_5)
print(f"ln(r_5) = {ln_r_5}")
# The contribution to the limit from this family.
lim_5 = 2 / ln_r_5
print(f"Limit contribution = 2 / {ln_r_5} = {lim_5}")

print("\n--- Final Calculation ---")
# The total limit is the sum of the two contributions.
total_lim = lim_13 + lim_5
print(f"Total limit L = {lim_13} + {lim_5} = {total_lim}")

# We need to find the integer part of 10^4 * L.
final_value = 10000 * total_lim
print(f"The final value is 10^4 * L = 10000 * {total_lim} = {final_value}")

# The integer part of the final value.
integer_part = int(final_value)
print(f"\nThe final equation to calculate is 10^4 * (2/ln((13+3*sqrt(17))/2) + 2/ln((5+sqrt(21))/2))")
print(f"The result is {final_value}, and its integer part is {integer_part}.")
<<<20637>>>