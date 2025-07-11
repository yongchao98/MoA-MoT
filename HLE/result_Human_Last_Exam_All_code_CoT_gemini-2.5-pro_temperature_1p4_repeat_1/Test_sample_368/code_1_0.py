import math

# This script calculates the probability that the given series converges.

# Define the possible values for X, Y, and Z
x_values = list(range(-9, 0)) + list(range(1, 10))
y_values = range(10)
z_values = range(10)

favorable_count = 0
total_count = 0

# Iterate through all possible combinations of X, Y, and Z
for x in x_values:
    for y in y_values:
        for z in z_values:
            total_count += 1
            
            # Let u be the term used to calculate the common ratio r.
            # u = (100*X + 10*Y + Z)/100 + Z/10 = X + Y/10 + 11*Z/100
            u = x + y / 10.0 + (11 * z) / 100.0
            
            # The common ratio of the geometric series
            r = 20 * u**2 + 24 * u
            
            # The series converges if the absolute value of the common ratio is less than 1.
            if -1 < r < 1:
                favorable_count += 1

# Calculate the probability as the ratio of favorable cases to total cases.
if total_count > 0:
    # Simplify the fraction using the greatest common divisor (GCD).
    common_divisor = math.gcd(favorable_count, total_count)
    numerator = favorable_count // common_divisor
    denominator = total_count // common_divisor
    
    # The final equation for probability is P = Favorable / Total.
    # The output shows the numbers in this equation.
    print(f"The number of favorable cases for convergence is: {favorable_count}")
    print(f"The total number of possible cases is: {total_count}")
    print(f"The probability is the fraction: {favorable_count} / {total_count}")
    print(f"The simplified probability is: {numerator} / {denominator}")
else:
    print("No valid cases were found to calculate the probability.")
