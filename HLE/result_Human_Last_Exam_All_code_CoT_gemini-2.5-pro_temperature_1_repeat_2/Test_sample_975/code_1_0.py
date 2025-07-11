import sympy as sp

# Define symbolic variables
M0, Rp, R, r, theta = sp.symbols('M_0 R_p R r theta')
ir, itheta = sp.symbols('i_r i_theta', real=True)

# --- Region 1: Inside the shield (0 < r < R_p) ---

# Coefficient for the H field inside the shield
coeff_H1 = M0 * (2*Rp**3 + R**3) / (3*R**3)

# Vector components for H1
H1_r_expr = -coeff_H1 * sp.cos(theta)
H1_theta_expr = coeff_H1 * sp.sin(theta)

# --- Region 2: Between shield and conductor (R_p < r < R) ---

# Radial component of H2
H2_r_expr = - (sp.S(2)/3) * M0 * ( (Rp/R)**3 - (Rp/r)**3 ) * sp.cos(theta)

# Theta component of H2
H2_theta_expr = (sp.S(1)/3) * M0 * ( 2*(Rp/R)**3 + (Rp/r)**3 ) * sp.sin(theta)


# --- Print the results ---
print("The correct answer is B.")
print("\n" + "="*50)
print("Derived Magnetic Field H(r, theta)")
print("="*50 + "\n")

print("In the region 0 < r < R_p (inside the shield):")
print("H = [ Hr ] * i_r   +   [ H_theta ] * i_theta\n")
print(f"Hr      = {sp.pretty(H1_r_expr, use_unicode=False)}")
print(f"H_theta = {sp.pretty(H1_theta_expr, use_unicode=False)}")
print("\nIn vector form:")
# Re-create the expression in the desired format from the option
H1_final_expr = f"H = M_0 * ((2*R_p**3 + R**3)/(3*R**3)) * (-cos(theta) * i_r + sin(theta) * i_theta)"
print(H1_final_expr)


print("\n" + "="*50 + "\n")

print("In the region R_p < r < R (between shield and conductor):")
print("H = [ Hr ] * i_r   +   [ H_theta ] * i_theta\n")
print(f"Hr      = {sp.pretty(H2_r_expr, use_unicode=False)}")
print(f"H_theta = {sp.pretty(H2_theta_expr, use_unicode=False)}")
print("\nIn vector form:")
# Re-create the expression in the desired format from the option
H2_final_expr_r = f"H_r = - (2*M_0/3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r"
H2_final_expr_theta = f"H_theta = (M_0/3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"
print(H2_final_expr_r)
print(H2_final_expr_theta)
