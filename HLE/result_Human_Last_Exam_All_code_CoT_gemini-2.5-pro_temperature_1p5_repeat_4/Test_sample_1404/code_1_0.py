import math
from fractions import Fraction

# Step 1: Calculate the number of arrangements for the large package.
# Number of ways to arrange 9 distinct red candles in 9 horizontal positions.
n_large_red_arrangements = math.factorial(9)
# Number of ways to arrange 9 distinct green candles in 9 vertical positions.
n_large_green_arrangements = math.factorial(9)
# Total arrangements for the large package is the product of the two.
N_large = n_large_red_arrangements * n_large_green_arrangements

# Step 2: Calculate the number of arrangements for the small package.
# We have 16 distinct candles (8 red, 8 green) to arrange in 16 horizontal positions.
# The total number of arrangements is 16!.
N_small = math.factorial(16)

# Step 3: Compare the number of arrangements.
# We check if N_small is 1260 times N_large by calculating the ratio.
ratio = N_small / N_large
target_ratio = 1260

# Step 4: Print the results and the final conclusion.
print("--- Large Package ---")
print(f"Ways to arrange 9 red candles horizontally: 9! = {n_large_red_arrangements}")
print(f"Ways to arrange 9 green candles vertically: 9! = {n_large_green_arrangements}")
print(f"Total arrangements (N_large) = 9! * 9! = {n_large_red_arrangements} * {n_large_green_arrangements} = {N_large}\n")

print("--- Small Package ---")
print(f"Ways to arrange 16 distinct candles horizontally (N_small): 16! = {N_small}\n")

print("--- Comparison ---")
print(f"The question is: Is N_small = {target_ratio} * N_large?")
print(f"Let's calculate the ratio: N_small / N_large")
print(f"Ratio = {N_small} / {N_large}")

# To avoid floating point inaccuracies for comparison, we can use fractions
is_true = (Fraction(N_small, N_large) == Fraction(target_ratio))

# Print the final result clearly
if is_true:
    print(f"\nThe ratio is exactly {target_ratio}.")
    print("Result: Yes, it is true that the number of possible arrangements for the small packages is 1260 times greater than for the large packages.")
else:
    # Represent the exact ratio as a simplified fraction
    exact_ratio = Fraction(N_small, N_large)
    print(f"\nThe calculated ratio is {ratio:.4f}...")
    print(f"The exact ratio as a fraction is {exact_ratio.numerator}/{exact_ratio.denominator}, which is not equal to {target_ratio}.")
    print("Result: No, it is not true that the number of possible arrangements for the small packages is 1260 times greater than for the large packages.")
