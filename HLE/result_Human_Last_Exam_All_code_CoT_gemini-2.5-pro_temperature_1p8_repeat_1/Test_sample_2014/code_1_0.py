import sympy

# --- Step 1: Define symbolic variables and parameters ---
# Let c be the chord length. We can set c=1 as the ratio is independent of it.
c = 1.0
# Separation s = c/2
s = c / 2
# Ride height h = c/2
h = c / 2
# G1, G2 are the circulations (proportional to lift L1, L2)
# K represents the lift of an isolated aerofoil
G1, G2, K = sympy.symbols('G1 G2 K')

# --- Step 2: Calculate induced vertical velocities (downwash) ---
# The vertical velocity 'w' (downwash positive) at point p(xp, yp) induced
# by a vortex of strength Gamma at v(xv, yv) is:
# w = (Gamma / (2*pi)) * (xp - xv) / ((xp - xv)^2 + (yp - yv)^2)

# Coordinates:
# V1(G1): (0, h), V2(G2): (s, h)
# V1_img(-G1): (0, -h), V2_img(-G2): (s, -h)
# Control Point 1 (CP1): (0, h), Control Point 2 (CP2): (s, h)

# --- Downwash at CP1 from all other vortices ---
# From V2
w_2_on_1 = (G2 / (2 * sympy.pi)) * (0 - s) / ((0 - s)**2 + (h - h)**2)
# From V1_img (image of itself)
w_1i_on_1 = (-G1 / (2 * sympy.pi)) * (0 - 0) / ((0 - 0)**2 + (h - (-h))**2)
# From V2_img
w_2i_on_1 = (-G2 / (2 * sympy.pi)) * (0 - s) / ((0 - s)**2 + (h - (-h))**2)
w1 = w_2_on_1 + w_1i_on_1 + w_2i_on_1

# --- Downwash at CP2 from all other vortices ---
# From V1
w_1_on_2 = (G1 / (2 * sympy.pi)) * (s - 0) / ((s - 0)**2 + (h - h)**2)
# From V1_img
w_1i_on_2 = (-G1 / (2 * sympy.pi)) * (s - 0) / ((s - 0)**2 + (h - (-h))**2)
# From V2_img (image of itself)
w_2i_on_2 = (-G2 / (2 * sympy.pi)) * (s - s) / ((s - s)**2 + (h - (-h))**2)
w2 = w_1_on_2 + w_1i_on_2 + w_2i_on_2

# --- Step 3: Formulate and solve the system of lift equations ---
# The fundamental lift equations are:
# G1 = K - pi*c*w1
# G2 = K - pi*c*w2
# (where K = pi*V*c*alpha_geom)
eq1 = sympy.Eq(G1, K - sympy.pi * c * w1)
eq2 = sympy.Eq(G2, K - sympy.pi * c * w2)

print("The system of equations for the circulations G1 and G2 is:")
# Use .simplify() to make the equations more readable
print(f"Equation 1: {sympy.simplify(eq1)}")
print(f"Equation 2: {sympy.simplify(eq2)}\n")

# --- Step 4: Solve for G1 and G2, then find the ratio ---
solution = sympy.solve([eq1, eq2], (G1, G2))
G1_sol = solution[G1]
G2_sol = solution[G2]
lift_ratio = sympy.simplify(G1_sol / G2_sol)

# Extract the coefficients for the final printout
g1_coeff = sympy.simplify(G1_sol / K)
g2_coeff = sympy.simplify(G2_sol / K)

print("Solving the system gives:")
print(f"Circulation of front aerofoil, G1 = ({g1_coeff}) * K")
print(f"Circulation of rear aerofoil,  G2 = ({g2_coeff}) * K\n")

print(f"The final lift ratio L1/L2 is calculated as G1/G2:")
print(f"L1/L2 = ({g1_coeff} * K) / ({g2_coeff} * K)")
print(f"L1/L2 = {lift_ratio}")
