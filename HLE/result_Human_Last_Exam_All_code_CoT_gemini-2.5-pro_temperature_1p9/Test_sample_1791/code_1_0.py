# Define the symbols used in the final expression
g = "g"
q_j = "q_j"
m = "m"
w_q = "w_q"
rho_q = "rho_q"
rho_nq = "rho_{-q}"

# Construct and print the equation for the effective interaction
# The equation is: H_eff = - (g^2 * q_j^2) / (2 * m * w_q^2) * rho_q * rho_{-q}

print("The effective electron-electron interaction for a single (q,j) mode is:")
print(f"H_eff = - ( {g}^2 * {q_j}^2 ) / ( 2 * {m} * {w_q}^2 ) * {rho_q} * {rho_nq}")

# The user asked to output each number/variable in the final equation.
# Here we print the components of the coefficient part of the equation.
print("\nComponents of the coefficient:")
print("Constant part: -1")
print(f"Coupling strength squared: {g}^2")
print(f"Momentum component squared: {q_j}^2")
print(f"Inverse of particle mass: 1/{m}")
print(f"Inverse of phonon frequency squared: 1/{w_q}^2")
print("Inverse of numerical factor: 1/2")
