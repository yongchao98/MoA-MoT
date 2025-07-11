import numpy as np
from scipy.special import lambertw

# We are solving the complex equation z*i = i^z for z in C.
#
# Let z = x + iy. The equation can be written as:
# -y + ix = exp((x+iy) * ln(i))
#
# The multi-valued logarithm of i is ln(i) = i * (pi/2 + 2*k*pi) for any integer k.
# Let C_k = pi/2 + 2*k*pi. The equation becomes:
# -y + ix = exp((x+iy) * i * C_k)
# -y + ix = exp(-y*C_k + i*x*C_k)
# -y + ix = exp(-y*C_k) * (cos(x*C_k) + i*sin(x*C_k))
#
# Equating the real and imaginary parts gives a system of two equations:
# 1) -y = exp(-y*C_k) * cos(x*C_k)
# 2)  x = exp(-y*C_k) * sin(x*C_k)

print("Solving the complex equation: z * i = i^z")
print("==========================================")

# Case 1: Real solutions (y=0)
# By substituting y=0 into the system, we can show that z=1 and z=-1 are solutions.
print("1. Real Solutions (z = x):")
real_solution_1 = 1
real_solution_2 = -1
print(f"Found two real solutions: z = {real_solution_1} and z = {real_solution_2}\n")

# Case 2: Pure imaginary solutions (x=0)
# By substituting x=0 into the system, we get:
# 1) -y = exp(-y*C_k)
# 2)  0 = 0 (This equation is satisfied for any y)
# We need to solve -y = exp(-y * (pi/2 + 2*k*pi)).
# Analysis shows that solutions for y only exist for k < 0.
# Let k = -j for j = 1, 2, 3, ...
# Let D_j = pi*(2j - 1/2). The equation becomes -y = exp(y*D_j).
# Substituting y = -a (for a>0), we get a = exp(-a*D_j), which leads to:
# a * D_j * exp(a * D_j) = D_j
# The solution is a*D_j = W(D_j), where W is the Lambert W function.
# So, a = W(D_j) / D_j, and the solutions are z = iy = -ia.

print("2. Pure Imaginary Solutions (z = iy):")
print("Solutions are of the form z = -i * W(D_j) / D_j, where D_j = pi*(2*j - 1/2) for j=1, 2, ...")
print("Here are the first 5 pure imaginary solutions:")

# Calculate and print the first 5 pure imaginary solutions
for j in range(1, 6):
    # D_j is the constant in the equation a * exp(a * D_j) = 1
    D_j = np.pi * (2 * j - 0.5)
    
    # For a positive real argument, the principal branch W_0 gives the real solution.
    W_Dj = lambertw(D_j).real
    
    # a is the positive real number such that z = -i*a
    a = W_Dj / D_j
    
    # The final solution for z
    z = -1j * a

    print(f"\n----- j = {j} -----")
    print(f"The equation for the solution's magnitude is a * exp(a * D_j) = 1")
    print(f"The constant in the equation is D_j = pi*(2*{j} - 0.5) = {D_j:.4f}")
    print(f"The Lambert W function argument is D_j = {D_j:.4f}")
    print(f"The value of the W function is W(D_j) = {W_Dj:.4f}")
    print(f"This gives a = W(D_j) / D_j = {a:.4f}")
    print(f"The final solution is z_{j} = {z.real:.4f}{z.imag:+.4f}i")

print("\n==========================================")
print("3. General Complex Solutions (x!=0 and y!=0):")
print("Other complex solutions exist but cannot be easily expressed in a closed form.")
print("They can be found by using numerical solvers on the system of equations for different integers k.")