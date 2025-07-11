import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

# Define the probability of A winning for a given r and F
def p_A_wins_r(r, F):
    if F == 0:
        return 0.5
    if r < F / 2:
        return 1.0
    
    # Check for domain of arccos
    arg = np.sqrt(2 * r * F - F**2) / r
    if arg > 1.0:
        arg = 1.0 # Should not happen based on analysis (r-F)^2 >= 0
    if arg < -1.0:
        arg = -1.0
        
    return np.arccos(arg) / np.pi

# Define the total probability P_A(F) by integrating over r
def P_A_total(F):
    if F == 0:
        return 0.5
    # The integral part
    integral, _ = quad(lambda r: p_A_wins_r(r, F), F/2, 1)
    # The first part of the probability
    prob = F / 2 + integral
    return prob

# Minimize P_A_total with respect to F
# We are looking for F in [0, 1]
res = minimize_scalar(P_A_total, bounds=(0, 1), method='bounded')

F_opt = res.x
P_min = res.fun

# Calculate the final result
final_value = int(np.floor(1 / P_min))

print(f"Optimal F: {F_opt}")
print(f"Minimized P(A wins): {P_min}")
print(f"1 / P(A wins): {1 / P_min}")
print(f"floor(1 / P(A wins)): {final_value}")

# To output the final equation clearly
# The final result is floor(1 / P_min)
# Equation: 1 / P_min
print(f"The equation we want to solve for is floor(1/P_min)")
print(f"The calculation leads to 1 / {P_min:.6f}")
final_val_text = 1/P_min
print(f"The result of this division is {final_val_text:.6f}")
final_answer = int(np.floor(final_val_text))
print(f"Taking the floor gives us {final_answer}")