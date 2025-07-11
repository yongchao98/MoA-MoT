import math

# Define the parameters from the problem description
A = 10**10
B = 1/100000 - 1
T_str = "10^20"

# The solution for X_0(t) is the constant A / (B + 1)
final_value = A / (B + 1)

# The final equation is X_0(T) = A / (B + 1)
# The prompt asks to output each number in the final equation.
# We will print the equation with the given numerical values.
print(f"X_0({T_str}) = {float(A)} / ({B} + 1) = {final_value:.0e}")