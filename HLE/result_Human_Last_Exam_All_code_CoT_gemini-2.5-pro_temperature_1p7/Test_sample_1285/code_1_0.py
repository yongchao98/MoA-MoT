# The problem asks for the maximum number of roots R_t can have in ]0, 1[.
# Based on the theory of Wronskians for functions of the form x^k * (1-x)^l,
# this maximum number is given by the formula C(t, 2) = t*(t-1)/2.

# (a) Provide the expression in terms of t.
expression_a = "t*(t-1)/2"

# (b) Calculate the value for t = 5.
t = 5
num1_b = t
num2_b = t - 1
result_b = num1_b * num2_b // 2

# Print the answers in the desired format.
print(f"(a) The expression for the maximum number of roots is: {expression_a}")

# As requested, output each number in the final equation for part (b).
print(f"(b) For t = {num1_b}, the calculation is:")
print(f"    {num1_b} * ({num1_b} - 1) / 2 = {num1_b} * {num2_b} / 2 = {result_b}")