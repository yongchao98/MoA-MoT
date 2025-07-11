import math

# Step 1: Define the given surface area.
area = 8.0

# Step 2: The theorem by Regina Rotman (2006) provides the upper bound for the
# length L of a closed geodesic on a 2-sphere with surface area A.
# The formula is L^2 <= 4 * pi * A, which simplifies to L <= 2 * sqrt(pi * A).

print("The upper bound for the length L of a closed geodesic is given by the formula:")
print("L <= 2 * sqrt(pi * A)")
print("")

# Step 3: Substitute the given area A = 8 into the formula.
print(f"Given A = {int(area)}, the equation becomes:")
print(f"L <= 2 * sqrt(pi * {int(area)})")
print("")

# Step 4: Simplify and calculate the final numerical value.
# The expression 2 * sqrt(8 * pi) can be simplified to 4 * sqrt(2 * pi).
print("The simplified form of the equation is:")
# Note: Here we show the mathematical simplification, not just the code variable
pi_symbol = "\u03C0" # Unicode for Greek letter pi
sqrt_symbol = "\u221A" # Unicode for square root symbol
print(f"L <= 4 * {sqrt_symbol}(2{pi_symbol})")
print("")

bound_value = 4 * math.sqrt(2 * math.pi)

print("The numerical value for this upper bound is:")
print(bound_value)
