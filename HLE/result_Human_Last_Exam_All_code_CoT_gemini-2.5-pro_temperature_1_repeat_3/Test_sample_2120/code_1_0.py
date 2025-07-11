import sys
from io import StringIO

# A brief explanation of the method within the code comments.
# The core of the problem is to find the average of the complex coordinates 'z'
# where the derivative of the Poynting-like field S(z) has poles.
# We deduce that these poles correspond to the roots of two specific polynomials
# derived from the equation for the B(z) field.

# Polynomial P(z) = 4z^4 - z^3 + z^2 + 1
# From the B(z) equation's right-hand side denominator.
# Coefficients: a=4, b=-1, c=1, d=0, e=1
P_coeffs = {'a': 4, 'b': -1}
num_roots_P = 4

# Polynomial Q(z) = z^4 + z^2 - z + 4
# From the denominator of the R(1/z) term when solving for B(z).
# Coefficients: a=1, b=0, c=1, d=-1, e=4
Q_coeffs = {'a': 1, 'b': 0}
num_roots_Q = 4

# According to Vieta's formulas, the sum of the roots of a polynomial
# a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0 is -a_{n-1}/a_n.

# Calculate sum of roots for P(z)
sum_roots_P = -P_coeffs['b'] / P_coeffs['a']

# Calculate sum of roots for Q(z)
sum_roots_Q = -Q_coeffs['b'] / Q_coeffs['a']

# The total set of coordinates is the union of the roots of P(z) and Q(z).
total_num_roots = num_roots_P + num_roots_Q
total_sum_of_roots = sum_roots_P + sum_roots_Q

# The average value is the total sum divided by the total number of roots.
average_z = total_sum_of_roots / total_num_roots

# Print the step-by-step calculation
print("The coordinates of interest are the roots of two polynomials, P(z) and Q(z).")
print("P(z) = 4z^4 - z^3 + z^2 + 1")
print("Q(z) = z^4 + z^2 - z + 4")
print(f"Total number of coordinates is the sum of the number of roots: {num_roots_P} + {num_roots_Q} = {total_num_roots}")
print(f"Sum of the roots of P(z) is -(-1)/4 = {sum_roots_P}")
print(f"Sum of the roots of Q(z) is -(0)/1 = {sum_roots_Q}")
print(f"The total sum of the coordinates is {sum_roots_P} + {sum_roots_Q} = {total_sum_of_roots}")
print(f"The final equation for the average value is ({sum_roots_P} + {sum_roots_Q}) / {total_num_roots}")
print(f"The average value of the complex coordinates z is {average_z}")

# Output the final answer in the specified format
# Create a string stream to capture the output for later formatting
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()
print(f"<<<{average_z}>>>")
sys.stdout = old_stdout
result = mystdout.getvalue()
