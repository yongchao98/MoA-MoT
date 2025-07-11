import sympy

# Define symbols
# G1, G2 are the circulations (Gamma) of the two aerofoils
# c is the chord length
# K is a constant representing the geometric angle of attack contribution (pi*c*U*alpha)
G1, G2, c, K = sympy.symbols('G1 G2 c K')

# Problem parameters
h = c/2  # Ride height
s = c/2  # Separation

# --- Vortex and Control Point Coordinates ---
# Aerofoil 1 (front)
x_v1, y_v1 = c/4, h
x_p1, y_p1 = 3*c/4, h

# Aerofoil 2 (rear)
x_v2, y_v2 = c + s + c/4, h
x_p2, y_p2 = c + s + 3*c/4, h

# Image vortices
x_v1_im, y_v1_im = x_v1, -h
x_v2_im, y_v2_im = x_v2, -h
G1_im_strength = -G1
G2_im_strength = -G2

# Define a function to calculate downwash 'w'
def get_downwash(vortex_G, vortex_pos, point_pos):
    """Calculates downwash at a point from a vortex."""
    vx, vy = vortex_pos
    px, py = point_pos
    r_sq = (px - vx)**2 + (py - vy)**2
    # The formula for downwash w from a CW vortex (positive lift)
    w = vortex_G * (px - vx) / (2 * sympy.pi * r_sq)
    return w

# --- Equation for Aerofoil 1 ---
# External downwash at P1 from G2, G1_im, G2_im
w_G2_on_P1 = get_downwash(G2, (x_v2, y_v2), (x_p1, y_p1))
w_G1im_on_P1 = get_downwash(G1_im_strength, (x_v1_im, y_v1_im), (x_p1, y_p1))
w_G2im_on_P1 = get_downwash(G2_im_strength, (x_v2_im, y_v2_im), (x_p1, y_p1))
w_ext_1 = w_G2_on_P1 + w_G1im_on_P1 + w_G2im_on_P1

# Equation 1: G1 = K - pi*c*w_ext_1
eq1 = sympy.Eq(G1, K - sympy.pi * c * w_ext_1)
eq1 = sympy.simplify(eq1)

# --- Equation for Aerofoil 2 ---
# External downwash at P2 from G1, G1_im, G2_im
w_G1_on_P2 = get_downwash(G1, (x_v1, y_v1), (x_p2, y_p2))
w_G1im_on_P2 = get_downwash(G1_im_strength, (x_v1_im, y_v1_im), (x_p2, y_p2))
w_G2im_on_P2 = get_downwash(G2_im_strength, (x_v2_im, y_v2_im), (x_p2, y_p2))
w_ext_2 = w_G1_on_P2 + w_G1im_on_P2 + w_G2im_on_P2

# Equation 2: G2 = K - pi*c*w_ext_2
eq2 = sympy.Eq(G2, K - sympy.pi * c * w_ext_2)
eq2 = sympy.simplify(eq2)

# Solve the system of equations for G1 and G2 in terms of K
solution = sympy.solve([eq1, eq2], (G1, G2))

# The lift ratio L1/L2 is equal to the circulation ratio G1/G2
lift_ratio = solution[G1] / solution[G2]

print("Equation for Aerofoil 1:")
# Sympy might rearrange the equation. We'll format it for clarity.
eq1_rhs = eq1.rhs.collect(G1).collect(G2)
eq1_final = sympy.Eq(eq1.lhs - eq1_rhs.coeff(G1)*G1, K + (eq1_rhs - K - eq1_rhs.coeff(G1)*G1)).simplify()
print(f"{eq1_final.lhs.coeff(G1)}*Γ1 + ({eq1_final.lhs.coeff(G2)})*Γ2 = K")


print("\nEquation for Aerofoil 2:")
eq2_rhs = eq2.rhs.collect(G1).collect(G2)
eq2_final = sympy.Eq(eq2.lhs - eq2_rhs.coeff(G2)*G2, K + (eq2_rhs - K - eq2_rhs.coeff(G2)*G2)).simplify()
#Need to handle coefficient being on RHS
term_g1 = eq2_final.rhs.coeff(G1)
term_k = eq2_final.rhs.coeff(K)
eq2_final_display = sympy.Eq(eq2_final.lhs - term_g1 * G1, term_k * K)

print(f"({eq2_final_display.lhs.coeff(G1)})*Γ1 + {eq2_final_display.lhs.coeff(G2)}*Γ2 = K")


print(f"\nSolving these two equations gives the lift ratio L1/L2 = Γ1/Γ2.")
print(f"L1/L2 = {lift_ratio}")