# The plan is to compute the sum based on the number of solutions for n odd/even.
# k_n is the number of possible values of m for a given n.
# For any n >= 1, m=0 and m=N_n=(3^n+(-1)^n)/2 are solutions. So k_n >= 2.
# For n odd, k_n = 2.
# For n even, there are two extra solutions, so k_n = 4.
#
# The total sum is Sum_{n=1 to inf} k_n / 4^n.
# We can write this as Sum_{n=1 to inf} 2/4^n + Sum_{n is even} 2/4^n.
#
# Part 1: Sum_{n=1 to inf} 2 * (1/4)^n
# This is a geometric series: 2 * ( (1/4)/(1-1/4) ) = 2 * ( (1/4)/(3/4) ) = 2 * (1/3) = 2/3.
#
# Part 2: Sum_{j=1 to inf} 2 * (1/4)^(2j) = Sum_{j=1 to inf} 2 * (1/16)^j
# This is a geometric series: 2 * ( (1/16)/(1-1/16) ) = 2 * ( (1/16)/(15/16) ) = 2 * (1/15) = 2/15.
#
# Total sum = 2/3 + 2/15 = 10/15 + 2/15 = 12/15 = 4/5.

from fractions import Fraction

# Calculate Part 1 of the sum
# Corresponds to 2 solutions for every n
part1_num = 2
part1_den = 3

# Calculate Part 2 of the sum
# Corresponds to 2 extra solutions for even n
part2_num = 2
part2_den = 15

total_sum_frac = Fraction(part1_num, part1_den) + Fraction(part2_num, part2_den)

final_num = total_sum_frac.numerator
final_den = total_sum_frac.denominator

# Output the equation and the final answer
print(f"The sum is given by the expression (2/3) + (2/15).")
print(f"({part1_num}/{part1_den}) + ({part2_num}/{part2_den}) = ({part1_num*5}/{part1_den*5}) + ({part2_num}/{part2_den}) = ({part1_num*5 + part2_num})/{part2_den} = {final_num}/{final_den}")
print(f"The final simplified answer is {total_sum_frac}.")
