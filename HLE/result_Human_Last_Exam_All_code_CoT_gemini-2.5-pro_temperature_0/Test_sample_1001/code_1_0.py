from fractions import Fraction

# The problem reduces to finding the sum of 1/2^(i+j) for three disjoint sets of pairs (i,j):
# 1. i = j, for k=1,2,3... (i=k, j=k)
# 2. j = 2i, for k=1,2,3... (i=k, j=2k)
# 3. i = 2j, for k=1,2,3... (i=2k, j=k)

# Case 1: i = j
# The sum is Sum_{k=1 to inf} 1/2^(2k) = Sum_{k=1 to inf} (1/4)^k
# This is a geometric series with a=1/4, r=1/4. Sum = a/(1-r) = (1/4)/(3/4) = 1/3.
sum1 = Fraction(1, 3)

# Case 2: j = 2i
# The sum is Sum_{k=1 to inf} 1/2^(k+2k) = Sum_{k=1 to inf} (1/8)^k
# This is a geometric series with a=1/8, r=1/8. Sum = a/(1-r) = (1/8)/(7/8) = 1/7.
sum2 = Fraction(1, 7)

# Case 3: i = 2j
# The sum is Sum_{k=1 to inf} 1/2^(2k+k) = Sum_{k=1 to inf} (1/8)^k
# This is the same sum as in Case 2.
sum3 = Fraction(1, 7)

# The total sum is the sum of the three cases.
total_sum = sum1 + sum2 + sum3

# The problem asks to output each number in the final equation.
print(f"The total sum is the sum of three series corresponding to the three cases:")
print(f"Case 1 (i=j): {sum1.numerator}/{sum1.denominator}")
print(f"Case 2 (j=2i): {sum2.numerator}/{sum2.denominator}")
print(f"Case 3 (i=2j): {sum3.numerator}/{sum3.denominator}")
print(f"Final Equation: {sum1.numerator}/{sum1.denominator} + {sum2.numerator}/{sum2.denominator} + {sum3.numerator}/{sum3.denominator} = {total_sum.numerator}/{total_sum.denominator}")
