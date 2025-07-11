import math

# --- Problem Parameters ---
# Half-width of the superconducting bar (m)
a = 0.01
# Half-thickness of the superconducting bar (m)
b = 0.0001
# Critical current density (A/m^2)
Jc = 3e10
# Applied magnetic field (A/m)
H = 20000

# --- Analytical Expression (as derived above) ---
# M(H) = -(a * Jc / 4) * tanh^2(π * H / (2 * b * Jc))

# --- Calculation ---
# 1. Calculate the effective critical sheet current, Kc = 2 * b * Jc
Kc = 2 * b * Jc

# 2. Calculate the argument of the hyperbolic tangent function
arg = (math.pi * H) / Kc

# 3. Calculate the magnetization prefactor
M_prefactor = -(a * Jc) / 4

# 4. Calculate the final magnetization using the formula
M = M_prefactor * (math.tanh(arg))**2

# --- Output ---
print("--- Analytical Expression ---")
print("M(H) = -(a * Jc / 4) * tanh^2(pi * H / (2 * b * Jc))")
print("\n--- Given Parameters ---")
print(f"a  = {a} m")
print(f"b  = {b} m")
print(f"Jc = {Jc:.1e} A/m^2")
print(f"H  = {H} A/m")

print("\n--- Equation with Numbers ---")
print(f"M({H}) = -({a} * {Jc:.1e} / 4) * tanh^2(pi * {H} / (2 * {b} * {Jc:.1e}))")

print("\n--- Step-by-step Calculation ---")
print(f"1. Prefactor -(a*Jc/4) = -({a} * {Jc:.1e} / 4) = {M_prefactor:.3e} A/m")
print(f"2. Argument (pi*H/(2*b*Jc)) = (pi * {H} / {Kc:.1e}) = {arg:.4f}")
print(f"3. M = {M_prefactor:.3e} * tanh^2({arg:.4f})")
print(f"4. M = {M_prefactor:.3e} * ({math.tanh(arg):.4f})^2")
print(f"5. M = {M_prefactor:.3e} * ({(math.tanh(arg))**2:.4f})")

print("\n--- Final Result ---")
print(f"The magnetization M is: {M:.4e} A/m")
<<<M(H) = -(a*Jc/4)*tanh^2(pi*H/(2*b*Jc))>>>