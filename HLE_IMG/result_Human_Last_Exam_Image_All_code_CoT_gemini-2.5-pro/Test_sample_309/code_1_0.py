import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def frobenius_number_two_vars(a1, a2):
    """Computes the Frobenius number for a set of two integers."""
    if gcd(a1, a2) != 1:
        return float('inf')
    return a1 * a2 - a1 - a2

# Based on the reasoning in the text:
# j=1
# m=51
# p=3
number_set = [51, 51 + 1, 3]
a1 = 3
a2 = 52 # 51 is a multiple of 3, so it's redundant.

# Compute the Frobenius number for {3, 52}
result = frobenius_number_two_vars(a1, a2)

print(f"The set of integers is {{{number_set[2]}, {number_set[0]}, {number_set[1]}}}.")
print(f"Since {number_set[0]} is a multiple of {number_set[2]}, the problem reduces to finding the Frobenius number of {{{number_set[2]}, {number_set[1]}}}.")
print(f"The formula for the Frobenius number of two variables {{a, b}} is a*b - a - b.")
print(f"Calculation: {a1} * {a2} - {a1} - {a2} = {a1*a2} - {a1} - {a2} = {result}")
print(f"The Frobenius number is {result}.")
