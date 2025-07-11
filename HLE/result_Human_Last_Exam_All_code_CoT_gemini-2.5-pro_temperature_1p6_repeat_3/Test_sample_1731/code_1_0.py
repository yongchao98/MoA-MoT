import math

# --- Introduction ---
print("This script presents the derivation of equilibrium mean energy and entropy for a photon gas.")
print("The foundation is Bose-Einstein statistics and the principle of maximum entropy, which is mathematically supported by Large Deviation Theory.")
print("-" * 70)

# --- Physical Constants (SI Units) ---
k_B = 1.380649e-23  # Boltzmann constant (J/K)
h = 6.62607015e-34   # Planck constant (J·s)
c = 2.99792458e8     # Speed of light (m/s)

# --- 1. Equilibrium Mean Energy <E> ---
print("1. EQUILIBRIUM MEAN ENERGY <E>\n")
print("The mean energy <E> of a photon gas in volume V at temperature T follows the Stefan-Boltzmann law for energy density.")
print("The symbolic formula is derived from the integral of energy over the Bose-Einstein distribution:")
print("    <E> = (8 * π⁵ * k_B⁴ / (15 * h³ * c³)) * V * T⁴\n")
print("This can be simplified to <E> = a * V * T⁴, where 'a' is the radiation constant.")

# Calculate the radiation constant 'a'
a_coeff_num = 8 * math.pi**5 * k_B**4
a_coeff_den = 15 * h**3 * c**3
a = a_coeff_num / a_coeff_den

print("To find the numerical value of the radiation constant 'a':")
print(f"    Numerator (8 * π⁵ * k_B⁴)   = {8:.1f} * {math.pi**5:.4f} * ({k_B:.4e})⁴ = {a_coeff_num:.4e}")
print(f"    Denominator (15 * h³ * c³) = {15:.1f} * ({h:.4e})³ * ({c:.4e})³ = {a_coeff_den:.4e}")
print(f"    a = Numerator / Denominator = {a:.4e} J·m⁻³·K⁻⁴\n")

print("The final equation for Mean Energy is:")
print(f"    <E> = {a:.4e} * V * T⁴  (in Joules)\n")
print("-" * 70)


# --- 2. Equilibrium Entropy S ---
print("2. EQUILIBRIUM ENTROPY S\n")
print("For a photon gas, the entropy S is related to the mean energy <E> by the simple thermodynamic relation:")
print("    S = (4/3) * <E> / T\n")

print("Substituting the formula for <E>, we get:")
print("    S = (4/3) * (a * V * T⁴) / T")
print("    S = (4/3 * a) * V * T³\n")
print("This can be simplified to S = b * V * T³, where 'b' is the entropy coefficient.")

# Calculate the entropy coefficient 'b'
b = (4.0 / 3.0) * a

print("To find the numerical value of the entropy coefficient 'b':")
print(f"    b = (4/3) * a = (4/3) * {a:.4e}")
print(f"    b = {b:.4e} J·m⁻³·K⁻³\n")

print("The final equation for Entropy is:")
print(f"    S = {b:.4e} * V * T³  (in Joules per Kelvin)\n")
