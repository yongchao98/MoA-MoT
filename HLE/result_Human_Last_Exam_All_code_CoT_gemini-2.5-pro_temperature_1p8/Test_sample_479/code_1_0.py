# The problem is a theoretical question from differential geometry.
# No computation is required to determine the answer, but we can state the known facts.
# Let g be the genus of the surface.

# Fact 1: A sphere has genus g=0.
# A sphere of radius R has constant mean curvature H = 1/R (or -1/R).
# This is non-vanishing. So, g=0 is a possible genus.
g_sphere = 0
print(f"A sphere is a surface of genus {g_sphere} with non-vanishing mean curvature.")

# Fact 2: A torus has genus g=1.
# A torus of revolution with major radius R and minor radius r has mean curvature
# that is always positive (and thus non-vanishing) if R > 2r.
# So, g=1 is a possible genus.
g_torus = 1
print(f"A torus can be constructed to have genus {g_torus} with non-vanishing mean curvature.")

# Fact 3: Higher genus surfaces (g > 1) cannot have non-vanishing mean curvature.
# This is a result of a theorem by Schoen and Yau.
# The theorem states that a compact, orientable, stable surface in a 3-manifold
# with non-negative scalar curvature has genus at most 1.
# - A surface with non-vanishing mean curvature (e.g., H > 0) is known to be stable.
# - R^3 has zero scalar curvature, so it is non-negative.
# - The surface in question is compact and orientable.
# Therefore, the theorem applies and the genus g must satisfy g <= 1.

max_g = 1
print(f"According to the Schoen-Yau theorem, the genus g must satisfy g <= 1.")
print(f"Since we have examples for g=0 and g=1, the maximal possible genus is {max_g}.")

final_answer = 'A'
# Final equation to show: maximal genus = 1
# Just for the output format requirement
equation_part_1 = "The maximal genus is"
equation_part_2 = "="
equation_part_3 = 1
print(f"{equation_part_1} {equation_part_2} {equation_part_3}")