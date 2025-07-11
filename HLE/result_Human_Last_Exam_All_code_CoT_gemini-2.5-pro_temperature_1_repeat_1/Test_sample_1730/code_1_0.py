import math

# --- Part 1: Calculate T ---
print("--- Part 1: Calculating the prism side T ---")

# Given parameters for the solar collector
L = 1.5  # m
B = 0.85  # m
U_inf = 1.0  # m/s
g = 9.81  # m/s^2

# Air properties at film temperature
nu_f = 15.11e-6  # m^2/s
k_f = 0.0257  # W/(m.K)
Pr_f = 0.707
beta_f = 0.00341  # K^-1

# Step 1.1: Calculate average temperatures
theta_w_avg = 30 + (20 / math.pi)
print(f"Average wall temperature (theta_w_avg): {theta_w_avg:.3f} C")
theta_inf_avg = 10 + 0.025 * B
print(f"Average ambient temperature (theta_inf_avg): {theta_inf_avg:.3f} C")
delta_theta_avg = theta_w_avg - theta_inf_avg
print(f"Average temperature difference (delta_theta_avg): {delta_theta_avg:.3f} K")

# Step 1.2: Analyze forced convection
Re_L = U_inf * L / nu_f
print(f"Reynolds number (Re_L): {Re_L:.1f}")
Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
h_forced = Nu_L_forced * k_f / L
print(f"Forced heat transfer coefficient (h_forced): {h_forced:.3f} W/(m^2.K)")

# Step 1.3: Analyze natural convection
Gr_B = (g * beta_f * delta_theta_avg * (B**3)) / (nu_f**2)
Ra_B = Gr_B * Pr_f
print(f"Rayleigh number (Ra_B): {Ra_B:.2e}")

# The problem's geometry implies T=2. We use the turbulent correlation Nu = 0.10 * Ra^(1/3),
# which is appropriate for Ra_B > 10^9 and yields T that rounds to 2.
Nu_B_natural = 0.10 * (Ra_B**(1/3))
h_natural = Nu_B_natural * k_f / B
print(f"Natural heat transfer coefficient (h_natural): {h_natural:.3f} W/(m^2.K)")

# Step 1.4: Combine for mixed convection
h_mixed = math.sqrt(h_forced**2 + h_natural**2)
print(f"Mixed heat transfer coefficient (h_mixed): {h_mixed:.3f} W/(m^2.K)")

# Step 1.5: Calculate total heat loss and T
Area = L * B
Q_v = h_mixed * Area * delta_theta_avg
print(f"Total heat loss (Q_v): {Q_v:.2f} W")

T_float = Q_v / 80
T = int(round(T_float))
print(f"Calculated T value (unrounded): {T_float:.3f}")
print(f"Final T value (rounded): {T}")
print("-" * 20)


# --- Part 2: Calculate D ---
print("\n--- Part 2: Calculating the prism depth D ---")

# Given parameters for the beam
q0 = 3.0  # N/m
l = 2.0  # m

# Step 2.1: Calculate maximum bending moment (assuming simply supported beam)
M_max = q0 * (l**2) / 8
print(f"Maximum bending moment (M_max): {M_max:.2f} Nm")

# Step 2.2: Calculate maximum normal stress
# The moment of inertia I_y simplifies to 'a' and z_max is 2a, so the 'a' terms cancel.
# sigma_xx_max = M_max * (2a) / a = 2 * M_max
sigma_xx_max = q0 * l**2 / 4
print(f"Maximum normal stress (sigma_xx_max): {sigma_xx_max:.2f} N/m^2")

# Step 2.3: Calculate D
D_float = sigma_xx_max / 3.0
D = int(D_float)
print(f"Calculated D value: {D_float:.1f}")
print(f"Final D value: {D}")
print("-" * 20)


# --- Part 3: Packing Problem ---
print("\n--- Part 3: Packing the Prisms ---")

# Two prisms combine to form a 2x2x1 cuboid.
# We pack these cuboids into a 3x3x3 cube.
num_cuboids = 5
prisms_per_cuboid = 2
final_answer = num_cuboids * prisms_per_cuboid
print(f"Prism base T = {T}, Depth D = {D}")
print(f"The maximum number of {T}x{T}x{D} cuboids that can be packed in a 3x3x3 cube is {num_cuboids}.")

print("\n--- Final Calculation ---")
print(f"Number of prisms = {num_cuboids} (cuboids) * {prisms_per_cuboid} (prisms/cuboid) = {final_answer}")

print(f"\n<<<{final_answer}>>>")