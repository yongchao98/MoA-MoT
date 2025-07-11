from fractions import Fraction

# Define the initial conditions
t0 = 0
x0 = -8
# Define the target time
t_target = 1

print("We are solving the initial value problem: x'(t) = (t-1) * x^2(t), with x(0) = -8.")
print("The equation is separable and its general solution is of the form: -1/x = t^2/2 - t + C\n")

# Step 1: Use the initial condition x(0) = -8 to find the constant C
print(f"Step 1: Find the constant C using the initial condition x({t0}) = {x0}.")
print("Substitute t = 0 and x = -8 into the general solution:")
# Equation: -1/x0 = t0^2/2 - t0 + C
C_lhs = Fraction(-1, x0)
C_rhs_term1 = Fraction(t0**2, 2)
C_rhs_term2 = -t0
C = C_lhs - C_rhs_term1 - C_rhs_term2

print(f"  -1/({x0}) = ({t0}^2)/2 - {t0} + C")
print(f"  {C_lhs} = {C_rhs_term1} - {t0} + C")
print(f"  From this, we find C = {C}\n")

# Step 2: Write the particular solution with the calculated C
print("Step 2: The particular solution is obtained by substituting C back into the general solution.")
print(f"  -1/x(t) = t^2/2 - t + {C}\n")

# Step 3: Calculate x(1) by substituting t=1
print(f"Step 3: Find x(1) by substituting t = {t_target} into the particular solution.")
# Equation: -1/x(1) = t_target^2/2 - t_target + C
rhs_val_at_t_target = Fraction(t_target**2, 2) - Fraction(t_target) + C
final_x_val = -1 / rhs_val_at_t_target

print(f"  -1/x({t_target}) = ({t_target}^2)/2 - {t_target} + {C}")
print(f"  -1/x({t_target}) = {Fraction(t_target**2, 2)} - {t_target} + {C}")
print(f"  -1/x({t_target}) = {rhs_val_at_t_target}")
print(f"Solving for x({t_target}):")
print(f"  x({t_target}) = {final_x_val}\n")

print(f"The final answer is x(1) = {final_x_val}.")
