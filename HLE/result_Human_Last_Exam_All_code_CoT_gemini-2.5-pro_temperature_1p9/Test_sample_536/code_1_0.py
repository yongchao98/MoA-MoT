import math

# Given values
alpha = 0.9375
beta = 0.9

# Coefficients of the quadratic equation x^2 + bx + c = 0
b = -2 * alpha * beta
c = alpha**2 + beta**2 - 1

# Calculate the discriminant
discriminant = b**2 - 4*c

print(f"The equation to solve is x^2 - 2*alpha*beta*x + (alpha^2 + beta^2 - 1) = 0")
print(f"Substituting alpha = {alpha} and beta = {beta}:")
print(f"x^2 + ({b:.4f})x + ({c:.8f}) = 0")
print(f"Discriminant = {discriminant:.8f}")

# Solve for the two roots
if discriminant < 0:
    print("No real solutions exist.")
else:
    x1 = (-b + math.sqrt(discriminant)) / 2
    x2 = (-b - math.sqrt(discriminant)) / 2
    print(f"The two possible solutions for x are: {x1:.6f} and {x2:.6f}")
    
    # In the context of these problems, a single value is expected.
    # Without further information to distinguish between the roots, and since
    # both roots are valid cosines, there might be ambiguity.
    # However, since the problem asks for a single limit, and no further constraints are given,
    # let's output the smaller of the two positive roots as it's a common convention in some contexts.
    # We choose one to provide a final single answer as requested.
    final_answer = x2
    
    print("\nThe final equation is:")
    # We chose the smaller root, x2.
    print(f"({final_answer:.6f})^2 - 2 * {alpha} * {beta} * {final_answer:.6f} + ({alpha})^2 + ({beta})^2 - 1 = 0")
    print(f"{final_answer**2:.6f} - {2 * alpha * beta * final_answer:.6f} + {alpha**2 + beta**2 - 1:.6f} = 0")
    
    print(f"\n<<< {final_answer} >>>")