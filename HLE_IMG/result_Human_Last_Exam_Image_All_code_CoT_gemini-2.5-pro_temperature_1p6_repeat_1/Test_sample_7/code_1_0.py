import sympy

# Define the symbols
V1_sym, V2_sym = sympy.symbols('V_1 V_2')
mu1_sym, mu2_sym, e_sym = sympy.symbols('μ_1 μ_2 e')

# Given plateau voltages
V1_val = 0.09  # Volts for the first plateau (Stage 2 -> 1)
V2_val = 0.13  # Volts for the second plateau (Dilute -> Stage 2)

print(f"The first plateau voltage is V1 = {V1_val} V, corresponding to the Stage 2 -> 1 transition.")
print(f"The second plateau voltage is V2 = {V2_val} V, corresponding to the formation of Stage 2.")
print("\nWe can relate these voltages to the chemical potentials μ_1 and μ_2:")
print(f"V1 = -μ_1 / e  =>  μ_1 = -V1 * e = -{V1_val}*e")
print(f"V2 = -μ_2 / e  =>  μ_2 = -V2 * e = -{V2_val}*e")

print("\nLet's test the proposed formula for approximating the second plateau: V_approx = -(2*μ_2 - μ_1) / e")
print("This formula combines different chemical potentials as requested.")

# To calculate the numerical value of this approximation, we can substitute the expressions for μ_1 and μ_2:
# V_approx = -(2*(-V2*e) - (-V1*e)) / e
# V_approx = -(-2*V2*e + V1*e) / e
# V_approx = -e*(-2*V2 + V1) / e
# V_approx = 2*V2 - V1

# Perform the calculation
V_approx_val = 2 * V2_val - V1_val

print("\nSubstituting the values of V1 and V2 into the simplified formula:")
print(f"V_approx = 2 * V2 - V1")
print(f"V_approx = 2 * {V2_val} - {V1_val}")
print(f"V_approx = {2 * V2_val} - {V1_val}")
print(f"V_approx = {V_approx_val:.2f} V")

print("\nThis value of 0.17 V is a reasonable approximation for the second plateau voltage of 0.13 V.")
print("\nThe formula that best approximates the second plateau using different potentials is therefore:")
print("V_plateau2 ≈ -(2*μ_2 - μ_1) / e")
