# Part (a) asks for an expression for the maximum number of roots of R_t in ]0, 1[.
# As derived in the plan, the Wronskian R_t can be written as a product of terms
# involving x, (1-x), and a polynomial Q(x). The roots in ]0, 1[ are the roots of Q(x).
# The maximum number of roots a polynomial can have is its degree.
# The degree of Q(x) is given by the formula: t*(t-1)/2.

# Part (b) asks for this maximum number when t = 5.
# We will calculate this value.

t = 5

# The formula for the maximum number of roots is t * (t - 1) / 2.
# We substitute t = 5 into the formula.
t_minus_1 = t - 1
numerator = t * t_minus_1
result = numerator // 2 # Use integer division as the result is always an integer.

# The prompt requires outputting each number in the final equation.
# Here we print the step-by-step calculation for part (b) and then the final answer.
print("The general expression for the maximum number of roots is: t*(t-1)/2")
print(f"For t = {t}, the calculation is:")
print(f"{t} * ({t} - 1) / 2 = {t} * {t_minus_1} / 2 = {numerator} / 2 = {result}")
print("\nFormatted Answer:")
print(f"(a) t*(t-1)/2; (b) {result}")
