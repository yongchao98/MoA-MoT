import math

# Define the symbols and terms for the equation
d_alpha = "\u2202\u2090" # ∂ₐ
sqrt_g = "\u221A(-g)" # √(-g)
P_sup_alpha_munu = "P\u1d43\u207B\u1d58\u1d5b" # Pᵃᵐᵘⁿᵘ
P_sub_mualphabeta = "P\u1d64\u2090\u1d66" # Pₘₐᵦ
Q_sup_nu_alphabeta = "Q\u1d65\u207B\u1d43\u1d66" # Qₙᵘᵃᵝ
Q_sup_alphabeta_mu = "Q\u1d43\u1d66\u1d64" # Qᵃᵝₘ
P_sub_alphabetanu = "P\u2090\u1d66\u1d65" # Pₐᵦₙ
Q_scalar = "Q"
g_munu = "g\u1d64\u1d65" # gₘₙ
T_munu = "T\u1d64\u1d65" # Tₘₙ

# Assemble the left-hand side (LHS) and right-hand side (RHS) of the equation
term1 = f"-2 / {sqrt_g} * {d_alpha}({sqrt_g} * {P_sup_alpha_munu})"
term2 = f"- {P_sub_mualphabeta} * {Q_sup_nu_alphabeta}"
term3 = f"+ 2 * {Q_sup_alphabeta_mu} * {P_sub_alphabetanu}"
term4 = f"+ 1/2 * {Q_scalar} * {g_munu}"

lhs = f"{term1} {term2} {term3} {term4}"
rhs = f"8\u03c0G/c\u2074 * {T_munu}"

# Print the final equation
print("The derived field equation is:")
print(f"{lhs} = {rhs}")

print("\nWhere the numbers and factors in the equation are:")
print(f"Factor of the first term on LHS: -2")
print(f"Factor of the second term on LHS: -1")
print(f"Factor of the third term on LHS: +2")
print(f"Factor of the fourth term on LHS: +1/2")
print(f"Proportionality constant on RHS: 8\u03c0G/c\u2074")