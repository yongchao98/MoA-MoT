import numpy as np

def z1(t):
    """
    Calculates the value of the trajectory z1(t).
    z1(t) = cos(t)
    """
    return np.cos(t)

def y(t):
    """
    Calculates the value of the entanglement echo y(t).
    y(t) = sin(2t) / (cos(2t))^(3/2)
    """
    return np.sin(2 * t) / (np.cos(2 * t))**(3/2)

# Set the time value
t_val = np.pi / 8

# Calculate z1(t) and y(t) at t = pi/8
z1_val = z1(t_val)
y_val = y(t_val)

# Calculate the final result
result = (z1_val / y_val)**2

# Print the results step-by-step
print(f"The calculation is for (z1(t) / y(t))^2 at t = pi/8.")
print(f"1. The value of t is: {t_val:.6f}")
print(f"2. The value of z1(t) = cos(t) is: {z1_val:.6f}")
print(f"3. The value of y(t) = sin(2t) / (cos(2t))^(3/2) is: {y_val:.6f}")
print(f"4. The ratio z1(t) / y(t) is: {z1_val/y_val:.6f}")
print(f"5. The final result (z1(t) / y(t))^2 is: {result:.6f}")

# For verification, the exact symbolic result is (1 + sqrt(2))/4
exact_result = (1 + np.sqrt(2)) / 4
print(f"\nThe exact symbolic result is (1 + sqrt(2))/4, which evaluates to: {exact_result:.6f}")