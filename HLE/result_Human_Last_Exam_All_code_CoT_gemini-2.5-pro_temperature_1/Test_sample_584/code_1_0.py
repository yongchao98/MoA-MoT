import math

# Given parameters for the S-curve
r = 3
b = 9
false_positive_rate = 0.05

# We are solving for the similarity threshold 's' in the equation:
# P(s) = 1 - (1 - s^r)^b
# where P(s) is the false positive rate.

# Rearranging the equation to solve for s:
# s = (1 - (1 - false_positive_rate)^(1/b))^(1/r)

# Perform the calculation
s_value = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

# --- Output the explanation and result ---

print("To find the threshold 's' where the false positive rate is 0.05, we solve the S-curve equation.")
print("The equation is: P(s) = 1 - (1 - s^r)^b")
print(f"With r={r}, b={b}, and P(s)={false_positive_rate}, we have:")
print(f"{false_positive_rate} = 1 - (1 - s^{r})^{b}\n")

print("Solving for 's' step-by-step:")
# Step 1
print("1. Rearrange the equation to isolate the term containing 's':")
val1 = 1 - false_positive_rate
print(f"   (1 - s^{r})^{b} = 1 - {false_positive_rate} = {val1}")

# Step 2
print("\n2. Take the {b}-th root of both sides:")
val2 = val1**(1/b)
print(f"   1 - s^{r} = ({val1})^(1/{b})")
print(f"   1 - s^{r} = {val2:.5f}")

# Step 3
print(f"\n3. Isolate s^{r}:")
val3 = 1 - val2
print(f"   s^{r} = 1 - {val2:.5f}")
print(f"   s^{r} = {val3:.5f}")

# Step 4
print(f"\n4. Take the {r}-th root to find 's':")
print(f"   s = ({val3:.5f})^(1/{r})")
print(f"   s = {s_value:.5f}")

# Final Answer
print("\n----------------------------------------------------")
print(f"The calculated threshold 's' is approximately {s_value:.5f}.")
print(f"Rounded to three decimal points, the threshold should be: {s_value:.3f}")
print("----------------------------------------------------")