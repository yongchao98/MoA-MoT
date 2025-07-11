# The exponent in the condition is 3/8
exponent_lambda = "3/8"

# We use the p=2 moment in Chebyshev's inequality.
# The measure of the set |X| is bounded by lambda^(-2) * integral.
# The exponent of N from this term is -2 * (3/8) = -3/4.
exponent_chebyshev = "-3/4"
numerator_chebyshev = -3
denominator_chebyshev = 4

# The integral term is bounded by the number of discretization points K.
# The number of points K needed is proportional to N^D, where D is the degree
# of the trigonometric polynomial in t. The frequencies are n^2, so the
# maximum frequency is N^2.
# So, the integral contributes a factor of N^2.
exponent_integral = "2"
integral_val = 2

# The final exponent alpha is the sum of these two exponents.
# alpha = 2 + (-3/4) = 2 - 3/4 = 5/4.
final_numerator = integral_val * denominator_chebyshev + numerator_chebyshev
final_denominator = denominator_chebyshev

print("The calculation for the exponent alpha is as follows:")
print(f"The threshold is of the order N^({exponent_lambda}).")
print(f"Using Chebyshev's inequality with p=2, we get a factor of N^(-2 * {exponent_lambda}) = N^({exponent_chebyshev}).")
print(f"The integral term is bounded by O(N^2), contributing a factor of N^({exponent_integral}).")
print(f"Thus, the final exponent alpha is {exponent_integral} + ({exponent_chebyshev}).")
print(f"alpha = {integral_val} - {abs(numerator_chebyshev)}/{denominator_chebyshev} = {final_numerator}/{final_denominator}")

# The final answer
alpha = final_numerator / final_denominator
# print(f"The numerical value is {alpha}")