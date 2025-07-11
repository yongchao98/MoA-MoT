import math

# This script calculates the number of iterations required for the
# Optimized Schwarz method to converge for the 1D wave equation.

# Problem parameters
# c: Propagation speed of the wave.
# L: Total length of the physical domain [0, L].
# T: The final time up to which the solution must be correct.
# The subdomain boundaries 'a' and 'b' are not needed for this calculation.

# --- You can change these example values ---
c = 1.0
L = 10.0
T = 25.0
# -----------------------------------------

# The theory for the Optimized Schwarz method states that the number of iterations (K)
# needed for convergence up to time T is the maximum number of times a wave, traveling
# at speed c, must reflect off the physical boundaries (x=0, x=L) to propagate
# information throughout the domain. This is equivalent to the number of times the
# wave traverses the full length of the domain, L.

# The total distance a wave can travel in time T is c * T.
# The number of full traversals of length L is floor( (c*T) / L ).
# The iteration counter starts at 0 and is incremented for each full update.
# If K traversals are needed, K updates are required, and the final counter
# value will be K.

num_iterations = math.floor((c * T) / L)

print("--- Calculation of Schwarz Iterations for the 1D Wave Equation ---")
print(f"Given parameters:")
print(f"  Propagation speed (c) = {c}")
print(f"  Domain Length (L)     = {L}")
print(f"  Final Time (T)        = {T}")
print("\nThe required number of iterations K is calculated as:")
print(f"  K = floor(c * T / L)")
print(f"  K = floor({c} * {T} / {L})")
print(f"  K = floor({c * T} / {L})")
print(f"  K = floor({(c * T) / L})")
print(f"  K = {num_iterations}")
<<<2>>>