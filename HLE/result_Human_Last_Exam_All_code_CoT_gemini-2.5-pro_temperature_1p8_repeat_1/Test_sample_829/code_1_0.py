import math

def F(u, u_bar):
    """Calculates F(u, u_bar)"""
    return u * (1 - u)**2 * math.exp(-u_bar)

def F11(u, u_bar):
    """Calculates the second derivative F_11"""
    return (-4 + 6 * u) * math.exp(-u_bar)
    
def P(u):
    """Calculates the polynomial P(u)"""
    return -12*u**3 + 24*u**2 - 16*u + 4

# Values that maximize the expression
u_x = 0.0
u_x_plus_1 = 1.0
u_bar_x = 0.0
u_bar_x_plus_1 = 0.0

# Calculate each component of the expression at the point (x)
P_at_u_x = P(u_x)
exp_term_at_x = math.exp(-2 * u_bar_x)
F11_at_x = F11(u_x, u_bar_x)

# Calculate the difference u(x+1) - u(x)
u_diff = u_x_plus_1 - u_x

# Calculate the flux F at points x and x+1
F_at_x = F(u_x, u_bar_x)
F_at_x_plus_1 = F(u_x_plus_1, u_bar_x_plus_1)
F_diff = F_at_x_plus_1 - F_at_x

# The expression to maximize has two main parts
# Part 1: P(u(x)) * exp(-2ū(x)) * (u(x+1) - u(x))
part1 = P_at_u_x * exp_term_at_x * u_diff

# Part 2: F₁₁(u(x),ū(x)) * (F(u(x+1),ū(x+1)) - F(u(x),ū(x)))
part2 = F11_at_x * F_diff

# The total expression is the sum of the two parts
max_value = part1 + part2

# Print the components of the final calculation
print(f"The calculation is performed at u(x) = {u_x}, u(x+1) = {u_x_plus_1}, u_bar(x) = {u_bar_x}, u_bar(x+1) = {u_bar_x_plus_1}")
print("\nFinal Equation: [P(u(x)) * exp(-2*u_bar(x))] * [u(x+1) - u(x)] + [F_11(u(x), u_bar(x))] * [F(u(x+1), u_bar(x+1)) - F(u(x), u_bar(x))]")
print("\nComponent values:")
print(f"P(u(x))              = {P_at_u_x}")
print(f"exp(-2*u_bar(x))     = {exp_term_at_x}")
print(f"u(x+1) - u(x)        = {u_diff}")
print(f"F_11(u(x), u_bar(x)) = {F11_at_x}")
print(f"F(u(x+1),u_bar(x+1)) = {F_at_x_plus_1}")
print(f"F(u(x), u_bar(x))    = {F_at_x}")
print(f"F_diff               = {F_diff}")

# Print the final assembly of the equation
print("\nAssembling the final equation with numbers:")
print(f"({P_at_u_x:.2f} * {exp_term_at_x:.2f}) * ({u_diff:.2f}) + ({F11_at_x:.2f}) * ({F_at_x_plus_1:.2f} - {F_at_x:.2f})")
print(f"= ({P_at_u_x * exp_term_at_x:.2f}) * ({u_diff:.2f}) + ({F11_at_x:.2f}) * ({F_diff:.2f})")
print(f"= {part1:.2f} + {part2:.2f}")
print(f"= {max_value:.2f}")

print("\nMaximum value:")
print(int(max_value))