import math

# This script calculates the minimal order of the Picard-Fuchs differential equation, u_r(n),
# for the real-analytic periods on the Hamiltonian with potential V(q) = 1/2 * (q^2 - q^n).
# The formula for this order is u_r(n) = 2 * floor((n-1)/2).
# We compute this for n in the range from 3 to 12.

print("Calculation of u_r(n) for n = 3 to 12 using the formula: u_r(n) = 2 * floor((n-1)/2)\n")

# This list will store the final sequence of results.
results = []

# Loop through each value of n from 3 to 12.
for n in range(3, 13):
    # The equation is u_r(n) = 2 * floor((n - 1) / 2).
    # We break it down to show each number in the calculation.
    numerator = n - 1
    value_inside_floor = numerator / 2
    floor_value = math.floor(value_inside_floor)
    result = 2 * floor_value

    # Store the calculated value.
    results.append(result)

    # Print the full calculation for the current value of n,
    # showing each number involved in the equation as requested.
    print(f"For n={n}: u_r({n}) = 2 * floor(({n}-1)/2) = 2 * floor({numerator}/2) = 2 * floor({value_inside_floor}) = 2 * {floor_value} = {result}")

# Print a separator for clarity before the final result.
print("\n" + "-"*70)

# Print the final complete set of values as a list.
print(f"The complete sequence of values {{u_r(3), u_r(4), ..., u_r(12)}} is:")
print(results)