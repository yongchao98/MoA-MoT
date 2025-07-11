import math

# Step 1: Calculate arrangements for the large package
n_large_factorial = math.factorial(9)
n_large = n_large_factorial * n_large_factorial

print("Step 1: Calculating arrangements for the large package.")
print("The 9 distinct red candles can be arranged in 9! ways.")
print("The 9 distinct green candles can be arranged in 9! ways.")
print(f"The number of arrangements for the large package is 9! * 9! = {n_large_factorial} * {n_large_factorial} = {n_large}\n")


# Step 2: Calculate arrangements for the small package
n_small = math.factorial(16)

print("Step 2: Calculating arrangements for the small package.")
print("There are 16 distinct candles to be arranged in 16 positions.")
print(f"The number of arrangements for the small package is 16! = {n_small}\n")


# Step 3: Compare the two quantities
ratio = n_small / n_large
target_ratio = 1260

print("Step 3: Comparing the two quantities.")
print(f"We need to check if N_small / N_large = {target_ratio}.")
print(f"The equation is: {n_small} / {n_large} = {target_ratio}?")
print(f"The actual ratio is: {ratio:.4f}\n")


# Step 4: Conclusion
is_true = (ratio == target_ratio)

print("Conclusion:")
if is_true:
    print(f"Yes, the statement is true. The number of arrangements for the small package is exactly {target_ratio} times greater than for the large package.")
else:
    print(f"No, the statement is false. The number of arrangements for the small packages is {ratio:.4f} times greater than for the large packages, which is not equal to {target_ratio}.")
