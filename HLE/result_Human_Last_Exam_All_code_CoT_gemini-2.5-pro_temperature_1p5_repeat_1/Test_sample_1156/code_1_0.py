import math

# Calculate the normalization constant C = 1 / (2 * ln(2))
C = 1 / (2 * math.log(2))

# The final equation for the normalized density rho(x)
# rho(x) = C / (sqrt(x) + x)
# We will print the numbers that appear in this equation.

# Numbers from the constant C: 1, 2, 2
# Numbers from the function body: 1 (implicit coefficient of sqrt(x)), 1 (implicit coefficient of x)
num1 = 1
num2 = 2
num3 = 2
num4 = 1
num5 = 1

print(f"The normalised density can be approximated by the function rho(x) = N / (sqrt(x) + x)")
print(f"The normalization constant N is 1/(2*ln(2)) which is approximately {C:.4f}.")
print("The final equation is of the form: N / (A*sqrt(x) + B*x)")
print(f"Here, N = 1/({num2}*ln({num3})), A = {num4}, B = {num5}")
