import math

# Step 1: Define parameters
# From the analysis of Bessel function intersections:
a = 5  # for n = 10000
lamb = 2  # for n = -2000

# Step 2: Determine N
# From the analysis of the quadratic equation for intersections near x=0:
# The condition for at most one intersection (zero real roots) is 9/4 - 40/n < 0
# which gives n < 160/9 = 17.77...
# For positive integers n, this means n can be 1, 2, ..., 17.
N = 17

# Step 3: Calculate y3(x0)
# From the solution of the fractional differential equation:
# y3(x0) = pi^(4.5) / 320
pi = math.pi
x0 = (pi / lamb)**lamb
# The calculation shows y3(x0) is pi**4.5 / 320
y3_x0 = (pi**4.5) / 320

# Step 4: Final calculation
# The expression to be calculated is (N + lambda) * (y3(x0))^(lambda/a)
final_value = (N + lamb) * (y3_x0)**(lamb / a)

print("This script calculates the solution to the multi-step physics problem.")
print("The key parameters were determined to be:")
print(f"a = {a}")
print(f"lambda = {lamb}")
print(f"N = {N}")
print("\nThe value of y3(x0) was calculated as:")
print(f"y3(x0) = pi^4.5 / 320 = {y3_x0}")
print("\nThe final expression is:")
# Print the equation with all the numbers.
# The user instruction is "Remember in the final code you still need to output each number in the final equation!"
print(f"({N} + {lamb}) * ({y3_x0})**({lamb}/{a})")
print("\nEvaluating this expression gives the final answer:")
print(f"({N+lamb}) * ({y3_x0})**({lamb/a}) = {final_value}")

print(f"\nFinal Answer: {final_value:.8f}")

# Wrapping the final numerical value in the requested format
print(f"\n<<<{final_value:.5f}>>>")