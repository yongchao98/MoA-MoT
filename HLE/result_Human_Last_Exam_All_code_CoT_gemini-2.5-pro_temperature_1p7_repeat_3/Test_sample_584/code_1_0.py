import math

# Given parameters for the S-curve
r = 3
b = 9
false_positive_rate = 0.05

# The equation is: P(s) = 1 - (1 - s^r)^b
# We need to find the threshold 's' for P(s) = 0.05
# So, we solve the equation: false_positive_rate = 1 - (1 - s^r)^b

# Derivation:
# 1 - (1 - s**r)**b = false_positive_rate
# (1 - s**r)**b = 1 - false_positive_rate
# 1 - s**r = (1 - false_positive_rate)**(1/b)
# s**r = 1 - (1 - false_positive_rate)**(1/b)
# s = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

threshold = (1 - (1 - false_positive_rate)**(1/b))**(1/r)

print("To find the threshold 's' for a false positive rate of 0.05, we solve the S-curve equation:")
print("P(s) = 1 - (1 - s^r)^b")
print("\nWith r=3, b=9, and P(s)=0.05, the equation to solve is:")
print(f"{false_positive_rate} = 1 - (1 - s**{r})**{b}")

print(f"\nThe calculated threshold 's' is: {threshold:.3f}")
print("\nAny similarity value below this threshold will have a false positive rate of less than 0.05.")

print("<<<0.178>>>")