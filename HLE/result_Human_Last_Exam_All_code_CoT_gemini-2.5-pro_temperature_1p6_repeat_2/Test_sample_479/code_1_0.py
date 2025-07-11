# This script formalizes the logical deduction based on theorems from differential geometry
# to find the maximal genus of the surface described in the problem.

# A theorem in differential geometry states that for a compact surface embedded in R^3,
# if its mean curvature never vanishes, its genus 'g' must be less than 2.
# This provides a strict upper bound.
upper_bound_for_genus = 2

# We must also show that a genus of 1 is achievable. A torus (genus 1) with major radius R
# and minor radius r has a mean curvature that is strictly positive everywhere if R > 2r.
# This confirms that the maximal genus is at least 1.

# The maximal genus must be an integer. Combining the two facts above, the genus 'g' must satisfy:
# g < 2 and g >= 1.
# The only integer that satisfies this is 1.

# We can express this conclusion as a simple equation.
# The highest possible integer value for the genus is the integer part of (upper_bound - epsilon),
# which is equivalent to upper_bound - 1 for an integer bound.
offset = 1
maximal_genus = upper_bound_for_genus - offset

# The problem requires printing the numbers in the final equation.
# Final equation: maximal_genus = upper_bound_for_genus - 1
print(f"A mathematical theorem imposes an upper bound on the genus (g). The bound is g < {upper_bound_for_genus}.")
print(f"Since genus 1 is known to be possible, the maximal genus is the largest integer satisfying the bound.")
print(f"Final Equation: {maximal_genus} = {upper_bound_for_genus} - {offset}")
print(f"The maximal genus is {maximal_genus}.")
