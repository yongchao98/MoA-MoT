import math

# This script derives and presents the analytical expression for the initial magnetization curve
# of a thin superconducting strip in a perpendicular magnetic field.

# --- Step 1 & 2: Define the physical model and parameters ---
# We consider a superconducting bar with a rectangular cross-section:
# Width = 2a (along x-axis, -a <= x <= a)
# Thickness = 2b (along y-axis, -b <= y <= b)
# Length = infinite (along z-axis)
# The problem assumes the thin strip limit, where b << a.
#
# An external magnetic field H is applied along the y-axis.
#
# The superconductor follows the critical-state model with a constant critical current density Jc.
# For a thin strip, we use the concept of a critical sheet current, I_c, which is the
# total current that can flow per unit width.
# I_c = Jc * 2b  (Amperes per meter)

# --- Step 3 & 4: Relate Applied Field H to Flux Penetration ---
# When the field H is applied, flux penetrates from the edges (x=±a) to a boundary (x=±c).
# In the penetrated region (c < |x| < a), a sheet current I(x) = -I_c * sgn(x) flows
# to screen the field. The current direction is chosen to create a magnetic moment opposing H.
#
# The magnetic field produced by this sheet current at the center of the strip (x=0) is H_ind.
# Using the Biot-Savart law for a strip, this field can be calculated as:
# H_ind(0) = - (I_c / pi) * ln(a/c)
#
# The screening condition requires the total field in the central region (|x|<c) to be zero:
# H_total = H + H_ind = 0
# H + (- (I_c / pi) * ln(a/c)) = 0
#
# This gives the relationship between the applied field H and the penetration boundary c:
# H = (I_c / pi) * ln(a/c)
# Rearranging for c/a, which we will need later:
# c/a = exp(-pi * H / I_c)

# --- Step 5: Calculate the Magnetization M ---
# Magnetization M is the magnetic moment per unit volume.
# The magnetic moment per unit length (m') is calculated by integrating the current distribution:
# m' = integral from -a to a of [x * I(x)] dx
# m' = integral from -a to -c of [x * I_c] dx + integral from c to a of [x * (-I_c)] dx
# m' = I_c * (c^2 - a^2)
#
# The volume per unit length is the cross-sectional area:
# Area = (2a) * (2b)
#
# Magnetization M = m' / Area = (I_c * (c^2 - a^2)) / (4ab)
# Substituting I_c = Jc * 2b:
# M = (Jc * 2b * (c^2 - a^2)) / (4ab)
# M = (Jc / (2a)) * (c^2 - a^2)
# We can rewrite this as:
# M = - (Jc * a / 2) * (1 - (c/a)^2)

# --- Step 6: Derive the Final Expression M(H) ---
# Now we substitute the expression for c/a from Step 4 into the expression for M from Step 5.
# c/a = exp(-pi * H / I_c)
# M = - (Jc * a / 2) * (1 - [exp(-pi * H / I_c)]^2)
# M = - (Jc * a / 2) * (1 - exp(-2 * pi * H / I_c))
#
# Finally, substitute I_c = Jc * 2b back into the expression:
# M = - (Jc * a / 2) * (1 - exp(-2 * pi * H / (Jc * 2b)))
# M = - (Jc * a / 2) * (1 - exp(-(pi * H) / (Jc * b)))

# --- Final Result ---
# The final analytical expression for the initial magnetization curve M(H) is:
print("The analytical expression for the initial magnetization M as a function of the applied field H is:")
print("")
print("  M(H) = M_sat * (1 - exp(-H / H_char))")
print("")
print("where the components of the equation are:")
print("  Jc: The critical current density of the superconductor.")
print("  a: Half of the strip's width.")
print("  b: Half of the strip's thickness.")
print("  pi: The mathematical constant " + str(math.pi))
print("")
print("And the composite terms are:")
# Printing each term as requested by the prompt
print("  Saturation Magnetization, M_sat = - (Jc * a / 2)")
print("  Characteristic Field, H_char = (Jc * b) / pi")
print("")
print("Combining these gives the full expression:")
print("  M(H) = - (Jc * a / 2) * (1 - exp(- (pi * H) / (Jc * b)))")
