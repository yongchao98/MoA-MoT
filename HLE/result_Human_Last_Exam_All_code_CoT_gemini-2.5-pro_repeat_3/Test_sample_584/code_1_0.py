import math

# Given parameters
r = 3
b = 9
false_positive_rate = 0.05

# Equation to solve for s: s = (1 - (1 - P)^(1/b))^(1/r)
# Where P is the desired probability (false_positive_rate)

print(f"Given r = {r} and b = {b}, we want to find the threshold 's' where the S-curve value is {false_positive_rate}.")
print("The equation to solve is: 1 - (1 - s^r)^b = P")
print("Rearranging for s, we get: s = (1 - (1 - P)^(1/b))^(1/r)\n")

print("Plugging in the values r=3, b=9, and P=0.05:")
print(f"s = (1 - (1 - {false_positive_rate})^(1/{b}))^(1/{r})\n")

# Step-by-step calculation
step1 = 1 - false_positive_rate
print(f"Step 1: 1 - {false_positive_rate} = {step1}")

step2 = step1**(1/b)
print(f"Step 2: {step1}^(1/{b}) = {step2}")

step3 = 1 - step2
print(f"Step 3: 1 - {step2} = {step3}")

s = step3**(1/r)
print(f"Step 4: {step3}^(1/{r}) = {s}\n")

# Final answer rounded to three decimal points
s_rounded = round(s, 3)
print(f"The threshold 's' should be approximately {s_rounded}.")

print(f"<<<{s_rounded}>>>")