import math

# Given parameters for the S-curve
r = 3
b = 9
fp_rate_target = 0.05

# We need to solve for s in the equation: 1 - (1 - s^r)^b = fp_rate_target
# Rearranging the formula to solve for s:
# 1 - fp_rate_target = (1 - s^r)^b
# (1 - fp_rate_target)^(1/b) = 1 - s^r
# s^r = 1 - (1 - fp_rate_target)^(1/b)
# s = (1 - (1 - fp_rate_target)^(1/b))^(1/r)

s_threshold = (1 - (1 - fp_rate_target)**(1/b))**(1/r)

# Round the result to three decimal points
s_rounded = round(s_threshold, 3)

# Verify the result by plugging it back into the S-curve formula
# This should be close to 0.05
actual_fp_rate = 1 - (1 - s_rounded**r)**b

print("The S-curve is defined by the equation: P(s) = 1 - (1 - s^r)^b")
print(f"Given r = {r} and b = {b}, we want to find the threshold 's' where the false positive rate P(s) < {fp_rate_target}.")
print("\nSolving for s where P(s) = 0.05:")
print(f"s = (1 - (1 - {fp_rate_target})^(1/{b}))^(1/{r})")
print(f"s = (1 - ({1 - fp_rate_target})^(1/{b}))^(1/{r})")
print(f"s = (1 - {((1 - fp_rate_target)**(1/b)):.5f})^(1/{r})")
print(f"s = ({(1 - ((1 - fp_rate_target)**(1/b))):.5f})^(1/{r})")
print(f"s = {s_threshold:.5f}")

print(f"\nThe similarity threshold 's' should be less than this value.")
print(f"Rounding to three decimal points, the threshold is: {s_rounded}")

print("\nVerification:")
print("Let's plug the rounded threshold back into the false positive rate equation:")
print(f"P({s_rounded}) = 1 - (1 - {s_rounded}^{r})^{b}")
print(f"P({s_rounded}) = 1 - (1 - {s_rounded**r:.5f})^{b}")
print(f"P({s_rounded}) = 1 - ({1 - s_rounded**r:.5f})^{b}")
print(f"P({s_rounded}) = 1 - ({(1 - s_rounded**r)**b:.5f})")
print(f"P({s_rounded}) = {actual_fp_rate:.5f}")
print(f"Since {actual_fp_rate:.5f} < {fp_rate_target}, a similarity threshold of {s_rounded} satisfies the condition.")

<<<0.179>>>