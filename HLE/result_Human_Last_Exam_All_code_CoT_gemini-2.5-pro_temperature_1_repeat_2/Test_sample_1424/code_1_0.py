import cmath

# --- Define Parameters ---
# omega: The frequency (energy) of the propagating particle.
omega = 2.5  # units of energy (e.g., eV)

# epsilon_k: A single-particle energy eigenvalue of the non-interacting system.
epsilon_k = 2.0  # units of energy (e.g., eV)

# delta: A small positive value to handle the pole.
# In a real calculation, this would be an infinitesimal limit.
delta = 0.1

# --- Calculation ---
# The formula for the bare Green's function is G_0 = 1 / (omega - epsilon_k + i*delta)
denominator = (omega - epsilon_k) + 1j * delta
G0 = 1 / denominator

# --- Output the result ---
print("The functional dependence of the bare Green's function G_0 on the single-particle energy eigenvalue ε_k is given by:")
print("G_0(k, ω) = 1 / (ω - ε_k + iδ)")
print("\nUsing the example values:")
print(f"  ω = {omega}")
print(f"  ε_k = {epsilon_k}")
print(f"  δ = {delta}")
print("\nThe final equation with these values is:")
print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{delta})")
print(f"G_0 = {G0}")