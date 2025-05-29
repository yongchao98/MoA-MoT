import math

successful_outcomes = 2665
total_outcomes = 10800

gcd = math.gcd(successful_outcomes, total_outcomes)
reduced_numerator = successful_outcomes // gcd
reduced_denominator = total_outcomes // gcd

print(f"{reduced_numerator}/{reduced_denominator}")