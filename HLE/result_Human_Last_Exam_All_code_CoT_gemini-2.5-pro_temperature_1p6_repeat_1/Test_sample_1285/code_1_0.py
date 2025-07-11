# Part (a): Determine the expression for the maximum number of roots.
t_symbol = 't'
a_expression = f"{t_symbol}*({t_symbol}-1)/2"

# Part (b): Calculate the maximum number of roots for t = 5.
t_value = 5
b_value = t_value * (t_value - 1) // 2

# To fulfill the requirement of showing the equation for the calculation,
# we format a string for part (b)'s answer.
b_expression_full = f"{t_value}*({t_value}-1)/2 = {b_value}"

# Print the final answer in the required format.
print(f"(a) The maximum number of roots is {a_expression}")
print(f"(b) For t = 5, this maximum number is {b_expression_full}")
