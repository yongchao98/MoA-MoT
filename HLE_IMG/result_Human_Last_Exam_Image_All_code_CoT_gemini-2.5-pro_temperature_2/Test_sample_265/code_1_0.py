import math

# --- Given Data ---
P = 90000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2
l_h = 50  # Length of horizontal welds W1 and W2 in mm
l_v = 100 # Length of vertical weld W3 in mm
ecc_dist = 150 # Distance from end of horizontal weld to force in mm

print("--- Step 1: Find the Centroid of the Weld Group ---")
# We set up a coordinate system with the origin at the bottom of the vertical weld (W3).
# Weld W1: from (0, 100) to (50, 100) -> Area A1=50, Centroid c1=(25, 100)
# Weld W2: from (0, 0) to (50, 0)     -> Area A2=50, Centroid c2=(25, 0)
# Weld W3: from (0, 0) to (0, 100)    -> Area A3=100, Centroid c3=(0, 50)
A1, c1x, c1y = l_h, l_h / 2, l_v
A2, c2x, c2y = l_h, l_h / 2, 0
A3, c3x, c3y = l_v, 0, l_v / 2
L = A1 + A2 + A3 # Total weld length

x_bar = (A1 * c1x + A2 * c2x + A3 * c3x) / L
y_bar = (A1 * c1y + A2 * c2y + A3 * c3y) / L

print(f"Total weld length, L = {A1} + {A2} + {A3} = {L} mm")
print(f"x_bar = (A1*c1x + A2*c2x + A3*c3x) / L = ({A1}*{c1x} + {A2}*{c2x} + {A3}*{c3x}) / {L} = {x_bar:.2f} mm")
print(f"y_bar = (A1*c1y + A2*c2y + A3*c3y) / L = ({A1}*{c1y} + {A2}*{c2y} + {A3}*{c3y}) / {L} = {y_bar:.2f} mm")
print(f"The centroid G of the weld group is at ({x_bar:.2f}, {y_bar:.2f}) mm.\n")

print("--- Step 2: Calculate Torque on the Weld Group ---")
# The force is applied at x = 50 + 150 = 200 mm from the coordinate system's y-axis.
x_force = l_h + ecc_dist
e = x_force - x_bar
T = P * e

print(f"Eccentricity, e = x_force - x_bar = {x_force} - {x_bar:.2f} = {e:.2f} mm")
print(f"Torque, T = P * e = {P} N * {e:.2f} mm = {T:.2f} N-mm\n")

print("--- Step 3: Calculate Unit Polar Moment of Inertia (Ju) ---")
# Ju = I_xx_u + I_yy_u, calculated about the centroid G.
# I_xx_u for weld 1 (top): A*dy^2 = 50 * (100-50)^2
I_xx_u1 = A1 * (c1y - y_bar)**2
# I_xx_u for weld 2 (bottom): A*dy^2 = 50 * (0-50)^2
I_xx_u2 = A2 * (c2y - y_bar)**2
# I_xx_u for weld 3 (vertical): l^3/12 + A*dy^2 = 100^3/12 + 100*(50-50)^2
I_xx_u3 = l_v**3 / 12 + A3 * (c3y - y_bar)**2
I_xx_u = I_xx_u1 + I_xx_u2 + I_xx_u3

# I_yy_u for weld 1 (top): l^3/12 + A*dx^2 = 50^3/12 + 50*(25-12.5)^2
I_yy_u1 = l_h**3 / 12 + A1 * (c1x - x_bar)**2
# I_yy_u for weld 2 (bottom): l^3/12 + A*dx^2 = 50^3/12 + 50*(25-12.5)^2
I_yy_u2 = l_h**3 / 12 + A2 * (c2x - x_bar)**2
# I_yy_u for weld 3 (vertical): A*dx^2 = 100 * (0-12.5)^2
I_yy_u3 = A3 * (c3x - x_bar)**2
I_yy_u = I_yy_u1 + I_yy_u2 + I_yy_u3

J_u = I_xx_u + I_yy_u

print(f"I_xx_u = ({I_xx_u1:.2f}) + ({I_xx_u2:.2f}) + ({I_xx_u3:.2f}) = {I_xx_u:.2f} mm^3")
print(f"I_yy_u = ({I_yy_u1:.2f}) + ({I_yy_u2:.2f}) + ({I_yy_u3:.2f}) = {I_yy_u:.2f} mm^3")
print(f"Unit polar moment of inertia, J_u = I_xx_u + I_yy_u = {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm^3\n")


print("--- Step 4: Determine Stress at Critical Point ---")
# The maximum stress occurs at the point farthest from the centroid where stress components add up.
# This is at the corners (50, 100) [Point A] and (50, 0) [Point B]. Let's analyze Point A.
px, py = l_h, l_v
# Radius vector from centroid G to point A
r_x = px - x_bar
r_y = py - y_bar
r = math.sqrt(r_x**2 + r_y**2)

print(f"Critical point A is at ({px}, {py}). It is the farthest from G and will have high stress.")
print(f"Vector from G to A: r_x = {px} - {x_bar:.2f} = {r_x:.2f} mm, r_y = {py} - {y_bar:.2f} = {r_y:.2f} mm.")

print("\n--- Step 5: Formulate Stress Equations in terms of Throat Thickness 't' ---")
# The force P creates a downward primary shear stress: τ'_y = P / (L * t)
# The torque T creates secondary shear stresses (τ''_x, τ''_y). For a clockwise torque:
# τ''_x = (T * r_y) / (J_u * t)
# τ''_y = -(T * r_x) / (J_u * t)
# Total stresses at Point A:
# τ_x = τ''_x
# τ_y = τ'_y + τ''_y = -(P/(L*t)) - (T*r_x)/(J_u*t) (both act downwards)

tau_primary_per_t = P / L
tau_secondary_x_per_t = (T * r_y) / J_u
tau_secondary_y_per_t = (T * r_x) / J_u

print(f"Primary shear stress: τ' = {P} / ({L} * t) = {tau_primary_per_t:.2f}/t")
print(f"Secondary shear stress (x-comp): τ''_x = ({T:.0f} * {r_y:.2f}) / ({J_u:.2f} * t) = {tau_secondary_x_per_t:.2f}/t")
print(f"Secondary shear stress (y-comp): τ''_y = -({T:.0f} * {r_x:.2f}) / ({J_u:.2f} * t) = -{tau_secondary_y_per_t:.2f}/t")

# Magnitudes of total stress components
tau_x_total_per_t = tau_secondary_x_per_t
tau_y_total_per_t = tau_primary_per_t + tau_secondary_y_per_t
print(f"\nTotal stress components at point A:")
print(f"τ_x = {tau_x_total_per_t:.2f} / t")
print(f"τ_y = - (τ' + τ''_y_mag) = -({tau_primary_per_t:.2f} + {tau_secondary_y_per_t:.2f}) / t = -{tau_y_total_per_t:.2f} / t\n")

print("--- Step 6: Solve for Throat Thickness (t) ---")
# τ_perm^2 = τ_x^2 + τ_y^2 = (1/t^2) * ( (τ_x_total_per_t)^2 + (τ_y_total_per_t)^2 )
# t^2 = ( (τ_x_total_per_t)^2 + (τ_y_total_per_t)^2 ) / τ_perm^2

t_squared = (tau_x_total_per_t**2 + tau_y_total_per_t**2) / tau_perm**2
t = math.sqrt(t_squared)

print(f"Using τ_perm^2 = τ_x^2 + τ_y^2")
print(f"{tau_perm}^2 = (1/t^2) * (({tau_x_total_per_t:.2f})^2 + ({tau_y_total_per_t:.2f})^2)")
print(f"{tau_perm**2} = (1/t^2) * ({tau_x_total_per_t**2:.2f} + {tau_y_total_per_t**2:.2f})")
print(f"t^2 = {t_squared:.4f} mm^2")
print(f"Required throat thickness, t = sqrt({t_squared:.4f}) = {t:.4f} mm\n")

print("--- Step 7: Solve for Weld Size (s) ---")
# For a fillet weld, t = 0.707 * s
s = t / 0.707

print(f"The relationship between leg size (s) and throat thickness (t) is t = 0.707 * s")
print(f"Weld leg size, s = t / 0.707 = {t:.4f} / 0.707 = {s:.4f} mm")

print(f"\nThe required size of the weld (leg) is {s:.2f} mm.")

<<<17.13>>>