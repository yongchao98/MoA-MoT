import math

# Calculate the initial position x(0)
# x(0) = 3 + (6*(3 - sqrt(3)))**(1/3) + (6*(3 + sqrt(3)))**(1/3)
# Let's simplify the terms inside the cube roots
sqrt3 = math.sqrt(3)
term1_cbrt_inner = 6 * (3 - sqrt3) # This is 18 - 6*sqrt(3)
term2_cbrt_inner = 6 * (3 + sqrt3) # This is 18 + 6*sqrt(3)

# Calculate the cube roots
term1 = term1_cbrt_inner**(1/3)
term2 = term2_cbrt_inner**(1/3)

# Calculate x(0)
x0 = 3 + term1 + term2

# The time t is given
t = 2 * sqrt3

# The trajectory is x(t) = x(0) + t^2 / 4
# Calculate the change in position
delta_x = t**2 / 4

# Calculate the final position x(t)
xt = x0 + delta_x

# The final result is the equation itself.
# x(t) = (3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)) + 3
# x(t) = 6 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
# We print the numbers in the final equation.

val_3 = 3
val_18 = 18
val_6_sqrt3 = 6 * sqrt3
final_t = t
val_t_squared_over_4 = delta_x
final_pos = xt

print("The initial position x(0) is calculated from the expression:")
print(f"x(0) = {val_3} + ({val_18:.4f} - {val_6_sqrt3:.4f})^(1/3) + ({val_18:.4f} + {val_6_sqrt3:.4f})^(1/3)")
print(f"x(0) = {val_3} + {term1_cbrt_inner:.4f}^(1/3) + {term2_cbrt_inner:.4f}^(1/3)")
print(f"x(0) ≈ {x0:.4f}\n")

print("The trajectory follows the equation x(t) = x(0) + t^2 / 4.")
print(f"At t = {final_t:.4f}, the displacement is t^2 / 4 = {val_t_squared_over_4:.4f}\n")

print("The final position x(t) is given by:")
print(f"x({final_t:.4f}) = x(0) + {val_t_squared_over_4:.4f}")
print(f"x({final_t:.4f}) ≈ {x0:.4f} + {val_t_squared_over_4:.4f} = {final_pos:.4f}")

# Final numerical answer
print("\nThe final numerical value is:")
print(f"{final_pos:.4f}")