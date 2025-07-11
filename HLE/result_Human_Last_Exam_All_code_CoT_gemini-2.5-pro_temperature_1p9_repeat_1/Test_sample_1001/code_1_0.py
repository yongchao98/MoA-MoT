from fractions import Fraction

# We identified three disjoint sets of pairs (i, j) that belong to S.
# 1. Pairs where i = j.
# 2. Pairs where j = 2i.
# 3. Pairs where i = 2j.
# We calculate the sum for each case.

# Case 1: i = j
# Sum = sum_{k=1 to inf} 1/2^(k+k) = sum_{k=1 to inf} (1/4)^k
# This is a geometric series with a=1/4, r=1/4. Sum = (1/4)/(1-1/4) = 1/3.
sum1 = Fraction(1, 3)

# Case 2: j = 2i
# Sum = sum_{i=1 to inf} 1/2^(i+2i) = sum_{i=1 to inf} (1/8)^i
# This is a geometric series with a=1/8, r=1/8. Sum = (1/8)/(1-1/8) = 1/7.
sum2 = Fraction(1, 7)

# Case 3: i = 2j
# Sum = sum_{j=1 to inf} 1/2^(2j+j) = sum_{j=1 to inf} (1/8)^j
# This is the same series as case 2. Sum = 1/7.
sum3 = Fraction(1, 7)

# The total sum is the sum of the sums from the three disjoint cases.
total_sum = sum1 + sum2 + sum3

print(f"The calculation is: {sum1.numerator}/{sum1.denominator} + {sum2.numerator}/{sum2.denominator} + {sum3.numerator}/{sum3.denominator} = {total_sum.numerator}/{total_sum.denominator}")
print(f"The numerical result is: {float(total_sum)}")

<<<13/21>>>