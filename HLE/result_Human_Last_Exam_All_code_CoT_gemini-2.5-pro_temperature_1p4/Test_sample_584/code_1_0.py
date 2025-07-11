import math

# Given parameters
r = 3
b = 9
# The target false positive rate
fp_rate = 0.05

# We are solving for 's' in the equation: P(s) = fp_rate
# P(s) = 1 - (1 - s^r)^b
#
# The equation becomes:
# fp_rate = 1 - (1 - s^r)^b
# 1 - fp_rate = (1 - s^r)^b
# (1 - fp_rate)^(1/b) = 1 - s^r
# s^r = 1 - (1 - fp_rate)^(1/b)
# s = (1 - (1 - fp_rate)^(1/b))^(1/r)

# --- Perform the calculation ---
val1 = 1 - fp_rate
val2 = val1**(1/b)
val3 = 1 - val2
threshold_s = val3**(1/r)

# --- Print the explanation and step-by-step solution ---
print("We want to find the similarity threshold 's' such that the probability P(s)")
print("of two items hashing to the same bucket is 0.05.")
print("The S-curve function is: P(s) = 1 - (1 - s^r)^b")
print(f"With r = {r} and b = {b}, the equation to solve is:")
print(f"1 - (1 - s^{r})^{b} = {fp_rate}")

print("\n--- Solving for s step-by-step ---")

# Print Step 1: Isolate the power term
print("\nStep 1: Rearrange the equation")
print(f"(1 - s^{r})^{b} = 1 - {fp_rate}")
print(f"(1 - s^{r})^{b} = {val1}")

# Print Step 2: Take the b-th root
print(f"\nStep 2: Take the {b}th root of both sides")
print(f"1 - s^{r} = ({val1})^(1/{b})")
print(f"1 - s^{r} = {val2:.6f}") # Show intermediate result

# Print Step 3: Isolate s^r
print(f"\nStep 3: Isolate s^{r}")
print(f"s^{r} = 1 - {val2:.6f}")
print(f"s^{r} = {val3:.6f}") # Show intermediate result

# Print Step 4: Take the r-th root
print(f"\nStep 4: Take the {r}rd root to solve for s")
print(f"s = ({val3:.6f})^(1/{r})")
print(f"s = {threshold_s:.6f}")

# Print final answer
print("\n--- Final Answer ---")
print(f"The threshold 's' should be {threshold_s:.3f} (rounded to three decimal places).")

print("\n<<<0.179>>>")