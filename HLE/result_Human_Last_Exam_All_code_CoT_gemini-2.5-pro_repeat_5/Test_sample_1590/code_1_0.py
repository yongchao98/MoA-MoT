import math

# This script derives the expression for the angle at which the rod begins to slide.
# The solution is symbolic as no numerical values are provided for the coefficient of friction.

# Plan:
# 1. Model the rod's interaction with the table corner using forces.
# 2. Resolve the force of gravity into components parallel and perpendicular to the rod.
# 3. Apply equilibrium conditions to find expressions for the Normal Force (N) and Friction Force (f).
# 4. Use the condition for impending motion (f = mu * N) to find the final equation and solve for theta.

print("Derivation of the angle (theta) at which the rod begins to slide:")
print("-" * 60)

# Step 1: Define forces in a coordinate system aligned with the rod.
print("1. The gravitational force (Mg) is resolved into two components based on the angle theta:")
print("   - Perpendicular to the rod: F_perp = M * g * cos(theta)")
print("   - Parallel to the rod (sliding force): F_para = M * g * sin(theta)")
print("")

# Step 2: Apply equilibrium conditions.
print("2. The normal force (N) from the table corner must balance the perpendicular component of gravity:")
print("   N = F_perp")
print("   N = M * g * cos(theta)")
print("")

print("3. To prevent sliding, the static friction force (f) must balance the parallel component of gravity:")
print("   f = F_para")
print("   f = M * g * sin(theta)")
print("")

# Step 3: Apply the condition for the onset of sliding.
print("4. Sliding begins when the required static friction force equals the maximum available static friction (f_max = mu * N):")
print("   M * g * sin(theta) = mu * N")
print("")

# Step 4: Substitute and solve for theta.
print("5. Substitute the expression for N from step 2 into the equation from step 4:")
print("   M * g * sin(theta) = mu * (M * g * cos(theta))")
print("")
print("6. The terms 'M' and 'g' (mass and gravity) cancel out, simplifying the equation:")
print("   sin(theta) = mu * cos(theta)")
print("")
print("7. The final equation is found by rearranging the terms to solve for theta:")
print("   sin(theta) / cos(theta) = mu")
print("   tan(theta) = mu")
print("")

# Final Answer
print("The final equation relating theta and the coefficient of friction mu is:")
print("tan(theta) = mu")
print("\nThe expression for the angle theta is therefore:")
print("theta = arctan(mu)")
print("\nComponents of the final equation:")
print("The term on the left side is: tan(theta)")
print("The term on the right side is: mu")