from fractions import Fraction

def geometric_series_sum(r_frac):
    """Calculates the sum of a geometric series from n=2 to infinity."""
    return (r_frac**2) / (1 - r_frac)

# Term for n=1
# For n=1, m can be 0 or 1. There are 2 pairs.
# Contribution is 2 * (1/4^1)
term_n1 = Fraction(2, 4)
print(f"Term for n=1: {term_n1}")

# Terms for n >= 2
# The sum is 1/2 * sum_{n=2 to inf} [2*(1/4)^n + (3/4)^n + (-1/4)^n]
r1 = Fraction(1, 4)
r2 = Fraction(3, 4)
r3 = Fraction(-1, 4)

sum1 = geometric_series_sum(r1)
sum2 = geometric_series_sum(r2)
sum3 = geometric_series_sum(r3)

sum_inf_part1 = 2 * sum1
sum_inf_part2 = sum2
sum_inf_part3 = sum3

print(f"Sum for n>=2 is composed of three parts:")
print(f"Part A (from 2*(1/4)^n): 2 * {sum1} = {sum_inf_part1}")
print(f"Part B (from (3/4)^n): {sum_inf_part2}")
print(f"Part C (from (-1/4)^n): {sum_inf_part3}")

sum_inf_total_inner = sum_inf_part1 + sum_inf_part2 + sum_inf_part3
sum_inf = Fraction(1, 2) * sum_inf_total_inner

print(f"The total infinite sum for n>=2 is: 1/2 * ({sum_inf_part1} + {sum_inf_part2} + {sum_inf_part3}) = 1/2 * {sum_inf_total_inner} = {sum_inf}")

# Total sum
total_sum = term_n1 + sum_inf
print(f"\nFinal Equation: {term_n1} + {sum_inf}")
print(f"The value of the sum is: {total_sum}")
