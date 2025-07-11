import math

# Given parameters
r = 3
b = 9
fp_rate = 0.05

# We need to find the threshold 's' by solving the equation for the boundary condition.
# The equation is: 1 - (1 - s^r)^b = fp_rate

print("The S-curve equation is: f(s) = 1 - (1 - s^r)^b")
print(f"We are given r = {r}, b = {b}, and we want the false positive rate to be < 0.05.")
print("To find the threshold, we solve for s when f(s) = 0.05:")
print(f"1 - (1 - s^{r})^{b} = {fp_rate}\n")

print("Step 1: Isolate the term with s.")
# 1 - fp_rate = (1 - s^r)^b
val1 = 1 - fp_rate
print(f"1 - {fp_rate} = (1 - s^{r})^{b}")
print(f"{val1} = (1 - s^{r})^{b}\n")

print(f"Step 2: Take the {b}-th root of both sides.")
# (1 - fp_rate)**(1/b) = 1 - s^r
val2 = val1**(1/b)
print(f"{val1}^(1/{b}) = 1 - s^{r}")
print(f"{val2:.6f} = 1 - s^{r}\n")

print(f"Step 3: Isolate s^{r}.")
# s^r = 1 - (1 - fp_rate)**(1/b)
val3 = 1 - val2
print(f"s^{r} = 1 - {val2:.6f}")
print(f"s^{r} = {val3:.6f}\n")

print(f"Step 4: Take the {r}-th root of both sides to solve for s.")
# s = (1 - (1 - fp_rate)**(1/b))**(1/r)
threshold = val3**(1/r)
print(f"s = ({val3:.6f})^(1/{r})")
print(f"s = {threshold:.6f}\n")

print("Final Answer:")
print(f"The threshold s, rounded to three decimal points, is {threshold:.3f}.")
print("For any similarity s < 0.178, the probability of being hashed to the same bucket is less than 0.05.")