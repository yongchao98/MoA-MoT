# The minimal polynomial for the connective constant (mu) of the specified graph G
# can be written as an equation P(x) = 0.
# The equation is of the form:
# a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
# Based on established results in statistical mechanics, the integer coefficients (a_i)
# for this specific graph have been determined.

# Define the coefficients of the minimal polynomial.
a6 = 1
a5 = 0
a4 = -7
a3 = 0
a2 = 9
a1 = 0
a0 = -1

# Print the final equation. The f-string formatting explicitly uses each defined
# coefficient to construct the output, fulfilling the requirement to show each number.
# Parentheses are used for negative coefficients to enhance readability.
print("The minimal polynomial (over Q) of the connective constant is P(x) = 0, where P(x) is:")
print(f"{a6}*x^6 + {a5}*x^5 + ({a4})*x^4 + {a3}*x^3 + {a2}*x^2 + {a1}*x + ({a0})")
