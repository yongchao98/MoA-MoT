import numpy as np

# Step 1: Analyze the system
# The system is X_{n+1} = A * X_n + b, where
# A = [[-1/2, -sqrt(3)/2], [sqrt(3)/2, -1/2]]
# A is a rotation matrix for an angle of 2*pi/3.
# A^3 is the identity matrix.
# This implies that the sequence X_n is periodic with period 3, i.e., X_{n+3} = X_n.
# We also need to show the expression we want to calculate is periodic.
# Let F_n = (x_n^1 - (3r/(2*pi))*sin(2*pi*n/3))^2 + (x_n^2)^2.
# Since sin(2*pi*(n+3)/3) = sin(2*pi*n/3) and X_{n+3}=X_n, F_{n+3} = F_n.

# Step 2: Simplify the target expression for n = 10^15
# We need to evaluate F_n for n = 10^15.
# We find the remainder of n divided by 3:
n = 10**15
n_mod_3 = n % 3
# Since 10 % 3 = 1, then 10^15 % 3 = 1^15 % 3 = 1.
# So, we need to calculate F_1.
# F_1 = (x_1^1 - (3r/(2*pi))*sin(2*pi/3))^2 + (x_1^2)^2
# F_1 = (x_1^1 - 3*sqrt(3)*r/(4*pi))^2 + (x_1^2)^2

# Let's relate F_1 to X_0.
# From the recurrence relation:
# x_1^1 = -1/2 * x_0^1 - sqrt(3)/2 * x_0^2 + 3*sqrt(3)*r/(4*pi)
# x_1^2 = sqrt(3)/2 * x_0^1 - 1/2 * x_0^2
# Substitute these into the expression for F_1:
# The first term becomes:
# ((-1/2*x_0^1 - sqrt(3)/2*x_0^2 + 3*sqrt(3)*r/(4*pi)) - 3*sqrt(3)*r/(4*pi))^2
# = (-1/2*x_0^1 - sqrt(3)/2*x_0^2)^2
# The second term is:
# (sqrt(3)/2*x_0^1 - 1/2*x_0^2)^2
# Adding them up and expanding:
# F_1 = (1/4*(x_0^1)^2 + 3/4*(x_0^2)^2 + sqrt(3)/2*x_0^1*x_0^2) + (3/4*(x_0^1)^2 + 1/4*(x_0^2)^2 - sqrt(3)/2*x_0^1*x_0^2)
# F_1 = (1/4 + 3/4)*(x_0^1)^2 + (3/4 + 1/4)*(x_0^2)^2 = (x_0^1)^2 + (x_0^2)^2.
# So, the value we need is ||X_0||^2.

# Step 3: Use boundary conditions to find X_0.
# First boundary condition: x_2025^2 = 10^20.
# Since 2025 is divisible by 3, X_2025 = X_0.
# So, x_0^2 = 10^20.

# Second boundary condition: -2/3*x_1^1 - 2/sqrt(3)*x_1^2 + sqrt(3)/(2*pi) = 0.
# Substitute the expressions for x_1^1 and x_1^2 in terms of x_0^1, x_0^2, and r:
# -2/3*(-1/2*x_0^1 - sqrt(3)/2*x_0^2 + 3*sqrt(3)*r/(4*pi)) - 2/sqrt(3)*(sqrt(3)/2*x_0^1 - 1/2*x_0^2) + sqrt(3)/(2*pi) = 0
# Simplifying this equation leads to:
# -2/3*x_0^1 + 2/sqrt(3)*x_0^2 - (sqrt(3)/(2*pi))*(r-1) = 0.

# Step 4: Determine r.
# For the problem to have a unique, well-defined numerical answer, the result cannot depend on the parameter r.
# The equation ` -2/3*x_0^1 + 2/sqrt(3)*x_0^2 - (sqrt(3)/(2*pi))*(r-1) = 0` must hold.
# The structure of the equation, particularly the `(r-1)` term, strongly suggests that `r=1`. With this assumption, the term involving `r` vanishes.
r = 1

# With r=1, the equation becomes:
# -2/3*x_0^1 + 2/sqrt(3)*x_0^2 = 0
# This simplifies to x_0^1 = sqrt(3)*x_0^2.

# Step 5: Calculate the final value.
x_0_2 = 10**20
x_0_1 = np.sqrt(3) * x_0_2

# The value to find is (x_0^1)^2 + (x_0^2)^2
result = x_0_1**2 + x_0_2**2

# We print the final equation for clarity.
x01_val_str = "(sqrt(3) * 10^20)"
x02_val_str = "(10^20)"
term1_val_str = "3 * 10^40"
term2_val_str = "10^40"
final_val_str = f"{result:.0e}"
# Reformat scientific notation
if 'e+' in final_val_str:
    base, exponent = final_val_str.split('e+')
    final_val_str = f"{base} * 10^{int(exponent)}"

print(f"The value is given by (x_0^1)^2 + (x_0^2)^2.")
print(f"Using the boundary conditions, we find x_0^1 = sqrt(3) * x_0^2 and x_0^2 = 10^20.")
print(f"So, the calculation is:")
print(f"{x01_val_str}^2 + {x02_val_str}^2 = {term1_val_str} + {term2_val_str} = {final_val_str}")
