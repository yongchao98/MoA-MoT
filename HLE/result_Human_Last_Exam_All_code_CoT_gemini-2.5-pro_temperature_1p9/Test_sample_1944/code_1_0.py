def get_js_number_str(n):
    """
    Creates a JSFuck-style string representation for a number.
    e.g., 2 -> '(! ![]+!![]) '
    """
    if n == 0:
        return '(+[])'
    if n == 1:
        return '(+! ![]) '
    if n > 1:
        return f"({''.join(['! ![]+'] * (n-1) ) }! ![]) "
    raise ValueError("Number must be non-negative.")

def get_js_string_concat(*args):
    """
    Creates a JS string concatenation expression.
    """
    return f"({'+'.join(args)})"

# In JavaScript, many of the expressions in the code evaluate to numbers.
# ! ![] is equivalent to the boolean `true`, and adding them together performs numeric addition.
# +!![] gives the number 1.
# +[] gives the number 0.

# 1. Let's define the numbers in our equation.
# The first number is 700. We can form the string "700" and then convert it to a number.
# The second number is 7.

# Representation of the number 7 in JSFuck style
js_seven = get_js_number_str(7) 
# Representation of the number 0 in JSFuck style
js_zero = get_js_number_str(0)

# 2. Construct the equation
# Form the string "700" by concatenating "7", "0", and "0"
# In JS, `'' + 7 + 0 + 0` results in the string `'700'`
js_str_700 = get_js_string_concat("''", js_seven, js_zero, js_zero)

# Convert the string '700' to a number in JS using the unary `+` operator
num1_str = f"+{js_str_700}"
num1_val = 700

# The second number is 7
num2_str = js_seven
num2_val = 7

# The final result
result = num1_val + num2_val

# 3. Print the final equation with each number.
print(f"The equation derived from the code is:")
print(f"{num1_val} + {num2_val} = {result}")
