import fractions

# We found the coefficients of f(x) to be:
a = fractions.Fraction(-3, 8)
b = fractions.Fraction(-5, 8)
c = fractions.Fraction(0, 1)

# The function is f(x) = x^3 + ax^2 + bx + c
# We need to compute f(3).
x = 3
result = x**3 + a * x**2 + b * x + c

# The problem asks to output each number in the final equation.
# f(3) = 3^3 + (-3/8)*3^2 + (-5/8)*3 + 0
val_x3 = x**3
val_ax2 = a * x**2
val_bx = b * x
val_c = c

print(f"The function is f(x) = x^3 + ({a})x^2 + ({b})x + {c}")
print(f"We want to compute f(3).")
print(f"f(3) = 3^3 + ({a})*3^2 + ({b})*3 + {c}")
print(f"f(3) = {val_x3} + ({val_ax2}) + ({val_bx}) + {c}")
print(f"f(3) = {val_x3} {val_ax2} {val_bx}")
# Combine the terms to show the calculation step by step
term_sum = val_ax2 + val_bx
print(f"f(3) = {val_x3} + ({term_sum})")
print(f"f(3) = {result}")

final_answer_str = str(result)
print(f"The exact value of f(3) is {final_answer_str}.")
print(f"<<<{final_answer_str}>>>")