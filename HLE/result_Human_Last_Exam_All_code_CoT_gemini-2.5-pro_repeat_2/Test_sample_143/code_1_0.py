import numpy as np

def g_double_prime(u, x1, x2):
    """
    Calculates the second derivative of the objective function g(u)
    where u = cos^2(theta).
    g(u) = (x1*sqrt(u) + x2*sqrt(1-u))^2
    The analytical formula for g''(u) is:
    g''(u) = -x1*x2*(1-u+u^2) / (2*(u-u^2)^(3/2))
    """
    if u <= 0 or u >= 1:
        return -np.inf if x1 > 0 and x2 > 0 else 0
    
    numerator = x1 * x2 * (1 - u + u**2)
    denominator = 2 * (u - u**2)**(1.5)
    
    if denominator == 0:
        return -np.inf if numerator > 0 else 0
        
    return -numerator / denominator

# --- Parameters for demonstration ---
# Let's choose a non-negative vector x in R^2
x = np.array([3.0, 4.0])
x1 = x[0]
x2 = x[1]

# Let's choose a value for u in (0, 1). u=0.5 corresponds to theta=pi/4.
u_val = 0.5

# --- Calculation ---
# Calculate the second derivative g''(u) at u_val
g_pp_val = g_double_prime(u_val, x1, x2)

# --- Output the results ---
print(f"Demonstration for d=2:")
print(f"Let x = ({x1}, {x2})")
print(f"We analyze the objective function g(u) = ({x1}*sqrt(u) + {x2}*sqrt(1-u))^2, where u = cos^2(theta).")
print(f"The second derivative g''(u) determines its concavity.")
print(f"The formula for the second derivative is: g''(u) = -({x1} * {x2} * (1 - u + u^2)) / (2 * (u - u^2)^1.5)")
print(f"\nEvaluating at u = {u_val}:")

# Manual calculation for clarity
numerator_val = x1 * x2 * (1 - u_val + u_val**2)
denominator_val = 2 * (u_val - u_val**2)**(1.5)

print(f"Numerator = {x1} * {x2} * (1 - {u_val} + {u_val**2}) = {numerator_val}")
print(f"Denominator = 2 * ({u_val} - {u_val**2})^1.5 = {denominator_val}")
print(f"g''({u_val}) = -{numerator_val} / {denominator_val} = {g_pp_val}")

print("\nSince g''(u) is negative, the objective function g(u) is concave.")
print("When all objective functions are concave functions of a single parameter, scalarization can find the entire Pareto front.")
print("This property holds for d=2 but fails for d>=3.")