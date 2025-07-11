from fractions import Fraction

# The optimal value for the parameter x0 was found to be 1/2.
x0 = Fraction(1, 2)

# Calculate the value of |b| for the optimal x0.
# Based on the formula: |b| = 4*x0 / (1+x0)^2
abs_b = 4 * x0 / (1 + x0)**2

# Calculate the value of |c| for the optimal x0.
# Based on the formula: |c| = (1 + 2*x0 - x0**2) / (1 + x0)**2
abs_c = (1 + 2*x0 - x0**2) / (1 + x0)**2

# The maximum value is the sum of |b| and |c|.
max_value = abs_b + abs_c

print("The final equation demonstrating the maximum value of |b| + |c| is:")
# The format requested is to output each number in the final equation.
print(f"{abs_b} + {abs_c} = {max_value}")
