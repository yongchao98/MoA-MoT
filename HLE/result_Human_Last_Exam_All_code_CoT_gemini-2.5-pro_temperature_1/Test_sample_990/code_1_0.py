import math

# The problem provides the height h in terms of g.
# h = 350g meters.
# Let's represent g symbolically for the explanation, as it will cancel out.
h_val = 350
h_expr = f"{h_val} * g"

# The derived formula for the time T to fall a distance h is:
# T = sqrt(14 * h / g)
print("The derived formula for the time T to fall a distance h is:")
print("T = sqrt(14 * h / g)")
print()

# Substitute the expression for h into the formula.
print(f"Substituting h = {h_expr} into the formula:")
print(f"T = sqrt(14 * ({h_expr}) / g)")
print()

# The 'g' terms cancel out, simplifying the equation.
print("The 'g' terms cancel out, simplifying the equation to:")
print(f"T = sqrt(14 * {h_val})")
print()

# Perform the multiplication inside the square root.
val_inside_sqrt = 14 * h_val
print(f"Calculating the value inside the square root: 14 * {h_val} = {val_inside_sqrt}")
print(f"T = sqrt({val_inside_sqrt})")
print()

# Calculate the final result by taking the square root.
final_time = math.sqrt(val_inside_sqrt)
print(f"The final time is the square root of {val_inside_sqrt}:")
print(f"T = {final_time}")
