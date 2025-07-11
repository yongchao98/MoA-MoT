from fractions import Fraction

# Step 1: Define initial conditions and target time
x0 = -8
t0 = 0
t1 = 1

# Step 2: Use the general solution (-1/x = t^2/2 - t + C) and initial conditions to find C.
# -1/x0 = t0^2/2 - t0 + C
# C = -1/x0 - (t0^2/2 - t0)
C = -Fraction(1, x0) - (Fraction(t0**2, 2) - t0)

print(f"The general solution is -1/x = t^2/2 - t + C.")
print(f"Using the initial condition x({t0}) = {x0}, we find C = {C}.")
print(f"The particular solution is: -1/x = t^2/2 - t + {C}\n")


# Step 3: Use the particular solution to find x(1).
# -1/x(1) = t1^2/2 - t1 + C
print(f"To find x({t1}), we solve the equation:")
val_t_sq_half = Fraction(t1**2, 2)
val_t = t1
print(f"-1/x({t1}) = ({t1})^2/2 - {t1} + {C}")
print(f"-1/x({t1}) = {val_t_sq_half} - {val_t} + {C}")

# Calculate the right-hand side
rhs = val_t_sq_half - val_t + C
print(f"-1/x({t1}) = {rhs}\n")

# Solve for x(1)
x1 = -1 / rhs
print(f"Therefore, x({t1}) = -1 / ({rhs}) = {x1}")

# Final numerical answer
# The result as a fraction is 8/3. As a float it is ~2.667
final_answer = float(x1)
print(f"\nThe final answer is {final_answer}")