from fractions import Fraction

# The problem is to solve x'(t) = (t-1)*x^2(t) with x(0) = -8 and find x(1).

# Initial condition values
t0 = 0
x0 = -8

# Find the constant C from the general solution -1/x = t^2/2 - t + C
# C = -1/x0 - (t0^2/2 - t0)
C = Fraction(-1, x0) - (Fraction(t0**2, 2) - Fraction(t0))

print("The general solution is of the form: -1/x(t) = t^2/2 - t + C")
print(f"Using the initial condition x({t0}) = {x0}, we find C:")
print(f"-1/({x0}) = ({t0})^2/2 - {t0} + C")
print(f"{Fraction(-1, x0)} = C")
print(f"So, the particular solution is: -1/x(t) = t^2/2 - t + {C}")
print("\nNow, we find x(1) by substituting t=1:")
print(f"-1/x(1) = (1)^2/2 - 1 + {C}")

# Calculate the right hand side of the equation for t=1
t = 1
rhs_val_1 = Fraction(t**2, 2)
rhs_val_2 = -t
rhs_val_3 = C
rhs_total = rhs_val_1 + rhs_val_2 + rhs_val_3

print(f"-1/x(1) = {rhs_val_1} - 1 + {rhs_val_3}")
print(f"-1/x(1) = {Fraction(4,8)} - {Fraction(8,8)} + {Fraction(1,8)}")
print(f"-1/x(1) = {rhs_total}")

# Solve for x(1)
x1 = -1 / rhs_total

print(f"Therefore, x(1) = 1 / ({-rhs_total})")
print(f"The final result is x(1) = {x1}")