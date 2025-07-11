import math

# Given parameters for the S-curve
r = 3
b = 9

# The target false positive rate
# This is the probability P(s) for a certain similarity threshold s
p = 0.05

# We need to solve for 's' in the equation: p = 1 - (1 - s^r)^b
print("The goal is to find the similarity threshold 's' given r, b, and a target probability p.")
print(f"Given values: r = {r}, b = {b}, p = {p}")
print("")
print("We start with the S-curve equation:")
print(f"P(s) = 1 - (1 - s^r)^b")
print("")
print("Substitute the given values to solve for s:")
print(f"{p} = 1 - (1 - s^{r})^{b}")
print("")
print("Step 1: Isolate the parenthetical term.")
print(f"(1 - s^{r})^{b} = 1 - {p}")
print(f"(1 - s^{r})^{b} = {1 - p}")
print("")
print("Step 2: Take the b-th root of both sides.")
print(f"1 - s^{r} = ({1-p})^(1/{b})")
# Calculate the intermediate value
val1 = (1-p)**(1/b)
print(f"1 - s^{r} = {val1}")
print("")
print("Step 3: Isolate the s^r term.")
print(f"s^{r} = 1 - {val1}")
# Calculate the next intermediate value
val2 = 1 - val1
print(f"s^{r} = {val2}")
print("")
print("Step 4: Take the r-th root of both sides to find s.")
print(f"s = ({val2})^(1/{r})")
# Calculate the final value of s
s = val2**(1/r)
print(f"s = {s}")
print("")
print("The similarity threshold 's' should be:")
print(f"{s:.3f}")
<<<0.179>>>