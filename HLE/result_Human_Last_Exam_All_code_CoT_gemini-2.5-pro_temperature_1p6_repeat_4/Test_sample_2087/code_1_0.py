# Define symbolic variables for the components of the equation
var_x = "x"
var_F_Xi_x = "F_{X_i}(x)"
var_mu_Xi = "\u03BC_{X_i}"  # Unicode for the Greek letter mu
var_I_Xi_x = "I_{X_i}(x)"

# The derived expression for the limiting CDF of X(t) is:
# (x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}

# The instruction is to output each "number" (or symbol) in the final equation.
# The following print statement constructs and displays the complete final equation,
# thereby showing all its constituent symbolic parts.

print("The final expression for lim_{t->\u221E} F_{X(t)}(x) is:")
print(f"({var_x} * {var_F_Xi_x} - {var_I_Xi_x}) / {var_mu_Xi}")