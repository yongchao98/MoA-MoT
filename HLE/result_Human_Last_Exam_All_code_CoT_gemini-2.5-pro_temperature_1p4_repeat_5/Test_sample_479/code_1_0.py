import math

# Step 1: Explain the problem for a torus (genus 1).
# The mean curvature H of a torus is non-vanishing if the term determining its sign,
# which is proportional to R + 2*r*cos(u), is never zero.
# This depends on the choice of the major radius R and minor radius r.

# Step 2: Choose parameters for a "thin" torus where R > 2r.
R = 3.0
r = 1.0

print(f"Let's test if a torus (genus 1) can have non-vanishing mean curvature.")
print(f"We choose a 'thin' torus with major radius R = {R} and minor radius r = {r}.")
print("The mean curvature H is zero if the following equation can be satisfied for some angle u:")

# Step 3: Output the equation with numbers, as requested.
print(f"{R} + 2 * {r} * cos(u) = 0")

# Step 4: Solve for cos(u) and check if a real solution exists.
# The equation can be rearranged to: cos(u) = -R / (2*r)
cos_u_for_zero_H = -R / (2 * r)
print(f"To solve this, we rearrange to find cos(u) = -{R} / (2 * {r}), which gives:")
print(f"cos(u) = {cos_u_for_zero_H}")

# Step 5: Analyze the result based on the possible range of cos(u).
# The function cos(u) can only take values between -1 and 1.
if abs(cos_u_for_zero_H) > 1:
    print(f"Since the required value {cos_u_for_zero_H} is outside the possible range [-1, 1] for cos(u), there is no real angle u that can solve the equation.")
    print("This means the mean curvature never vanishes for this specific torus.")
    print("Therefore, a surface of genus 1 with non-vanishing mean curvature is possible.")
else:
    # This branch is not executed for the chosen R and r.
    print("A solution for u exists, so the mean curvature vanishes at some points for this torus.")

# Step 6: State the concluding fact about higher genera.
print("\nA known, though deep, theorem in differential geometry states that any compact surface")
print("with genus g >= 2 must have at least one point where the mean curvature is zero.")

# Step 7: Final conclusion.
print("\nCombining these facts, since genus 1 is possible but any genus >= 2 is not, the maximal possible genus is 1.")