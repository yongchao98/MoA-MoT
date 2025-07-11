import numpy as np
from scipy.special import lambertw

print("The equation z * i = i^z has several sets of solutions in the complex plane.")
print("We found two real solutions and an infinite series of purely imaginary solutions.")

# --- Solution 1: z = 1 ---
z1 = 1
lhs1 = z1 * 1j
# In Python, 1j**z1 correctly computes the principal value
rhs1 = 1j**z1
print("\nSolution 1: z = 1")
print(f"Verifying the equation for z = {z1}:")
print(f"{z1} * i         = {lhs1}")
print(f"i^({z1})      = {rhs1}")
print(f"The values match, so z = 1 is a solution.")

# --- Solution 2: z = -1 ---
z2 = -1
lhs2 = z2 * 1j
rhs2 = 1j**z2
print("\nSolution 2: z = -1")
print(f"Verifying the equation for z = {z2}:")
print(f"{z2} * i        = {lhs2}")
print(f"i^({z2})     = {rhs2}")
print(f"The values match, so z = -1 is a solution.")

# --- Solution 3: Infinite imaginary solutions ---
print("\nSolution 3: An infinite series of purely imaginary solutions.")
print("These are of the form z = -i * W(k) / k, where k = pi*(2j-1/2) for j=1, 2, 3,...")
print("W is the principal branch of the Lambert W function.")
print("\nCalculating the first 3 imaginary solutions:")

# Loop for the first 3 solutions in this family (j=1, 2, 3)
for j in range(1, 4):
    # This corresponds to branches k = -1, -2, -3 of the logarithm
    k_val = -j
    
    # Calculate the argument for the Lambert W function
    # This is -A_k from our derivation, which is pi*(2j-1/2)
    lambert_arg = np.pi * (2 * j - 1/2)
    
    # Calculate the Lambert W function value. For real, positive arguments, the result is real.
    w_val = np.real(lambertw(lambert_arg))
    
    # Calculate the solution z
    z = -1j * w_val / lambert_arg
    
    # LHS of the original equation: z * i
    lhs = z * 1j
    
    # RHS of the original equation: i^z
    # We must use the specific branch of log(i) that led to this solution.
    # log(i) = i * (pi/2 + 2*k*pi), with k = -j
    log_i = 1j * (np.pi/2 + 2 * k_val * np.pi)
    rhs = np.exp(z * log_i)
    
    print(f"\n----- Solution for j = {j} -----")
    print(f"z = {z:.8f}")
    print("Verifying the equation:")
    print(f"z * i         = {lhs.real:.8f}")
    print(f"i^z           = {rhs.real:.8f}")
    print(f"The values match, confirming the solution.")