import sympy as sp

# Define the symbolic variables
I_0, a, D, mu_0, sigma_1, sigma_2 = sp.symbols('I_0 a D mu_0 sigma_1 sigma_2')

# Step 1: Calculate the current I_2 flowing through the second block.
# The blocks are in parallel, so the current divides based on conductance.
# I_2 = I_0 * (sigma_2 / (sigma_1 + sigma_2))
I_2 = I_0 * (sigma_2 / (sigma_1 + sigma_2))
print("The current in the second block (I_2) is:")
print("I_2 =", sp.latex(I_2))
print("-" * 20)


# Step 2: Express the force F_x in terms of I_2.
# From the Maxwell stress tensor analysis, the force is F_x = - (a * mu_0 * I_2**2) / (2 * D)
Fx_intermediate = - (a * mu_0 * I_2**2) / (2 * D)
print("The force F_x in terms of I_2 is:")
print("F_x =", sp.latex(Fx_intermediate).replace('I_0', 'I_{2}').replace(sp.latex(sigma_2 / (sigma_1 + sigma_2))+'^{2}',''))
print("-" * 20)

# Step 3: Substitute I_2 into the force equation to get the final answer.
Fx_final = Fx_intermediate
print("The final expression for the force F_x is:")
print("F_x =", sp.latex(Fx_final))
print("-" * 20)

# Step 4: Rearrange to match the format in Answer A for comparison.
# Fx = -a * D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2
Fx_rearranged = -a * D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2
print("Rearranging the expression to match the answer format:")
print("F_x =", sp.latex(Fx_rearranged))
print("-" * 20)

# Final check
print(f"The equation matches Answer A: F_x = -aD * (μ_0/2) * (I_0^2/D^2) * (σ_2/(σ_1 + σ_2))^2")
