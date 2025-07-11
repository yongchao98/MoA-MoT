# Define the numbers from the expression
base = 12
exponent = 4
# From the analytical solution, the final result is determined by this numerator
numerator = 5

# Calculate the value of base^exponent
term_A = base**exponent

# The problem is to compute (12)^4 * (Integral)^4.
# Our derivation shows Integral = (5^(1/4))/12.
# So, (Integral)^4 = ( (5^(1/4))/12 )^4 = 5 / (12^4).
# The full expression is (12^4) * (5 / 12^4), which equals 5.
# The following line prints this final equation with the computed numbers.
# It shows the equation in the form: A * (B / A) = B
print(f"{term_A} * ({numerator} / {term_A}) = {numerator}")