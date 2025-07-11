import math

# Given parameters for the S-curve
r = 3
b = 9
# The desired false positive rate
fp_rate = 0.05

# Explain the setup and the equation being solved.
print(f"We are given r = {r} and b = {b} for the S-curve P(s) = 1-(1-s^r)^b.")
print(f"We need to find the similarity threshold 's' where the false positive rate is {fp_rate}.")
print("\nThis means we must solve the following equation for 's':")
print(f"1 - (1 - s^{r})^{b} = {fp_rate}")

# Solving the equation for s:
# s = (1 - (1 - fp_rate)**(1/b))**(1/r)
s_value = (1 - (1 - fp_rate)**(1/b))**(1/r)

print("\nSolving for 's' step-by-step:")
print(f"(1 - s^{r})^{b} = 1 - {fp_rate}")
print(f"(1 - s^{r})^{b} = {1 - fp_rate}")
print(f"1 - s^{r} = ({1-fp_rate})**(1/{b})")
print(f"s^{r} = 1 - ({1-fp_rate})**(1/{b})")
print(f"s = (1 - ({1-fp_rate})**(1/{b}))**(1/{r})")

# Print the calculated value and the rounded result.
print(f"\nThe calculated value for s is: {s_value}")
s_rounded = round(s_value, 3)
print(f"The threshold 's' rounded to three decimal points is: {s_rounded}")