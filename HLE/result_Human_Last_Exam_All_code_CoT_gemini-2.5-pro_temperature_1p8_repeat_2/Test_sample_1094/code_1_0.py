# This script calculates and displays the normalized AC loss for an elliptic
# superconductor carrying a transport current, based on Norris's formula from
# the critical-state model.

# The AC loss per cycle per unit length (Q) for a normalized current i = Im/Ic < 1 is:
# Q = (μ₀ * Ic²) / π * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i²]

# We are asked to find the expression for 2*π*Q / (μ₀*Ic²).
# Let's perform the normalization:
# 1. Multiply Q by 2*π:
#    2*π*Q = 2*π * { (μ₀*Ic²) / π * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i²] }
#    The 'π' terms cancel out:
#    2*π*Q = 2 * μ₀*Ic² * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i²]
#
# 2. Divide by (μ₀*Ic²):
#    (2*π*Q) / (μ₀*Ic²) = { 2 * μ₀*Ic² * [...] } / (μ₀*Ic²)
#    The 'μ₀*Ic²' terms cancel out, leaving the final expression.

# The code below prints the final derived equation.
print("For an elliptic superconductor in the critical state carrying a transport current,")
print("the normalized loss per cycle per unit length Q is given by:")
print("")
# The following line prints the formula with all its numerical components, as requested.
print("2*π*Q / (μ₀*Ic²) = 2*((1-i)*ln(1-i) + (1+i)*ln(1+i) - i²)")
print("")
print("where:")
print("  i = Im/Ic is the transport current amplitude normalized to the critical current.")
print("  ln is the natural logarithm.")
