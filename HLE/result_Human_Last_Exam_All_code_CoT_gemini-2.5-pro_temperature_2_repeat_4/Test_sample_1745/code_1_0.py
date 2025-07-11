import math

# This script provides the final expression for the Electrical double-layer (EDL)
# potential distribution, ψ(y), under the specified conditions.

# --- Parameters in the Equation ---
# ψ(y): EDL potential as a function of position y
# z₁:   Base zeta potential at the bottom surface (y=0)
# β:    Slip length
# k:    Debye–Huckel parameter
# H:    Total height of the microchannel
# y:    Position coordinate across the channel height (from 0 to H)
# sinh: The hyperbolic sine function, defined as sinh(x) = (e^x - e^-x) / 2

print("The final expression for the Electrical double-layer potential distribution ψ(y) is:")
print("-" * 70)

# The full expression is derived from the linearized Poisson-Boltzmann equation
# with boundary conditions ψ(0) = z₁(1+βk) and ψ(H) = 0.
final_equation = "ψ(y) = z₁ * (1 + β * k) * sinh(k * (H - y)) / sinh(k * H)"

print(final_equation)
print("-" * 70)

print("\nWhere each term and number in the final equation represents:")
print(f"  ψ(y)  -> The potential at position y.")
print(f"  z₁    -> A constant representing the base zeta potential on the bottom surface.")
print(f"  1     -> The number one, originating from the slip-dependency term for zeta potential (1 + βk).")
print(f"  β     -> A constant representing the slip length.")
print(f"  k     -> The Debye-Hückel parameter, related to the ionic strength of the fluid.")
print(f"  H     -> The total height of the microchannel.")
print(f"  y     -> The variable for vertical position within the channel (from 0 to H).")
print(f"  sinh  -> The hyperbolic sine function.")