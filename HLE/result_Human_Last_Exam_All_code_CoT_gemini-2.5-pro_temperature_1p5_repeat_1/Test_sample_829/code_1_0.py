import numpy as np
from scipy.optimize import shgo

# Define the functions from the problem statement
def F_val(u, b):
    """Calculates F(u, bar_u) = u(1-u)^2 * exp(-bar_u)"""
    return u * (1 - u)**2 * np.exp(-b)

def P_val(u):
    """Calculates P(u) = 4 - 16u + 24u^2 - 12u^3"""
    return 4 - 16*u + 24*u**2 - 12*u**3

def F11_val(u, b):
    """Calculates F11(u, bar_u) = (-4+6u) * exp(-bar_u)"""
    return (-4 + 6*u) * np.exp(-b)

def expression_to_maximize(x):
    """
    Calculates the full expression E for a given state vector x.
    x = [u1, b1, u2, b2]
    u1=u(x), b1=bar_u(x), u2=u(x+1), b2=bar_u(x+1)
    """
    u1, b1, u2, b2 = x
    
    # The expression derived in the text:
    # P(u1)*exp(-2*b1)*(u2-u1) + F11(u1,b1)*(F(u2,b2) - F(u1,b1))
    
    val_P_u1 = P_val(u1)
    val_F11_u1_b1 = F11_val(u1, b1)
    val_F_u2_b2 = F_val(u2, b2)
    val_F_u1_b1 = F_val(u1, b1)

    term1 = val_P_u1 * np.exp(-2 * b1) * (u2 - u1)
    term2 = val_F11_u1_b1 * (val_F_u2_b2 - val_F_u1_b1)
    
    return term1 + term2

# Since SciPy optimizers are minimizers, we minimize the negative of our expression
def objective_function(x):
    return -expression_to_maximize(x)

# Define the bounds for the variables [u1, b1, u2, b2], which are all in [0, 1]
bounds = [(0, 1), (0, 1), (0, 1), (0, 1)]

# Use the shgo (simplicial homology global optimization) algorithm for global optimization
result = shgo(objective_function, bounds)

# The maximum value is the negative of the minimum found by the optimizer
max_value = -result.fun
# The point at which the maximum is achieved
max_point = result.x

u1_opt, b1_opt, u2_opt, b2_opt = max_point

# Output the results, including the values of each term at the maximum
print("The expression to maximize is E = (P(u1)*exp(-2*b1)*(u2-u1)) + (F11(u1,b1)*(F(u2,b2) - F(u1,b1)))")
print("\n--- Values at Maximum ---")
print(f"The maximum value of the expression is: {max_value:.4f}")
print(f"This maximum is achieved at the point (u1, b1, u2, b2) = ({u1_opt:.4f}, {b1_opt:.4f}, {u2_opt:.4f}, {b2_opt:.4f})")
print("\n--- Breakdown of the expression at the maximum point ---")
P_u1 = P_val(u1_opt)
F11_u1_b1 = F11_val(u1_opt, b1_opt)
F_u2_b2 = F_val(u2_opt, b2_opt)
F_u1_b1 = F_val(u1_opt, b1_opt)
term1 = P_u1 * np.exp(-2 * b1_opt) * (u2_opt - u1_opt)
term2 = F11_u1_b1 * (F_u2_b2 - F_u1_b1)

print(f"P(u1) = {P_u1:.4f}")
print(f"F11(u1, b1) = {F11_u1_b1:.4f}")
print(f"F(u2, b2) = {F_u2_b2:.4f}")
print(f"F(u1, b1) = {F_u1_b1:.4f}")
print(f"First term of E = {term1:.4f}")
print(f"Second term of E = {term2:.4f}")
print(f"Total E = First term + Second term = {term1 + term2:.4f}")

<<<4.0>>>