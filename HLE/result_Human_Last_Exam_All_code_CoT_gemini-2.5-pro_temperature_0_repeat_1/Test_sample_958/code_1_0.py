import sympy as sp

def solve_evanescent_energy():
    """
    This function symbolically derives the time-averaged stored energy per unit area
    for an evanescent wave created by a p-polarized light in TIR.
    """
    # 1. Define all symbols. We assume positive real values where appropriate.
    n, theta, omega, c, z = sp.symbols('n, theta, omega, c, z', real=True, positive=True)
    epsilon_0, mu_0 = sp.symbols('epsilon_0, mu_0', real=True, positive=True)
    # The problem gives the incident field in terms of the amplitude of its x-component.
    E_x0_i_amp = sp.Symbol('|E_{x0}^i|', real=True, positive=True)

    # 2. Define intermediate quantities based on wave optics for TIR.
    # The incident angle theta is greater than the critical angle, so n*sin(theta) > 1.
    sin_theta = sp.sin(theta)
    cos_theta = sp.cos(theta)
    
    # Wave vector component parallel to the interface
    k_x = n * (omega / c) * sin_theta
    
    # Decay constant perpendicular to the interface
    alpha = (omega / c) * sp.sqrt(n**2 * sin_theta**2 - 1)
    
    # Cosine of the complex transmission angle
    cos_theta_t = sp.I * sp.sqrt(n**2 * sin_theta**2 - 1)

    # 3. Calculate the squared amplitude of the transmitted electric field.
    # The transmission coefficient t_p relates the total E-field amplitudes.
    # t_p = E_transmitted / E_incident
    # The incident E-field amplitude is related to its x-component by E_incident = E_x0_i / cos(theta)
    t_p = (2 * n * cos_theta) / (1 * cos_theta + n * cos_theta_t)
    E_0_i_amp = E_x0_i_amp / cos_theta
    E_0_t_complex = t_p * E_0_i_amp
    
    # The squared magnitude of the transmitted E-field amplitude
    E_0_t_amp_sq = sp.simplify(sp.Abs(E_0_t_complex)**2)
    
    # Let's manually simplify the denominator for clarity, which sympy will also do.
    # denom = (n**2 - 1) * ((n**2 + 1) * sin_theta**2 - 1)
    # E_0_t_amp_sq = (4 * n**2 * E_x0_i_amp**2) / denom

    # 4. Define the complex field amplitudes in the cladding (z > 0).
    # The transmitted E-field has x and z components. The H-field has a y component.
    E_tx_complex_z = E_0_t_complex * cos_theta_t * sp.exp(-alpha * z)
    E_tz_complex_z = E_0_t_complex * (-n * sin_theta) * sp.exp(-alpha * z)
    
    Z_0 = sp.sqrt(mu_0 / epsilon_0) # Impedance of free space
    H_ty_complex_z = (1 / Z_0) * E_0_t_complex * sp.exp(-alpha * z)

    # 5. Formulate the time-averaged energy densities.
    # <u_E> = (1/4) * epsilon_0 * |E_t|^2 = (1/4) * epsilon_0 * (|E_tx|^2 + |E_tz|^2)
    # <u_H> = (1/4) * mu_0 * |H_t|^2
    u_E_z = sp.simplify((1/4) * epsilon_0 * (sp.Abs(E_tx_complex_z)**2 + sp.Abs(E_tz_complex_z)**2))
    u_H_z = sp.simplify((1/4) * mu_0 * (sp.Abs(H_ty_complex_z)**2))

    # 6. Integrate densities over z from 0 to infinity to get energy per unit area.
    # The integral of exp(-2*alpha*z) from 0 to oo is 1/(2*alpha).
    integral_factor = 1 / (2 * alpha)
    
    W_E = sp.simplify(u_E_z.subs(z, 0) * integral_factor)
    W_H = sp.simplify(u_H_z.subs(z, 0) * integral_factor)

    # 7. Print the final results.
    print("Derived Time-Averaged Stored Energy in the Electric Field (W_E):")
    sp.pretty_print(W_E)
    print("\n" + "="*50 + "\n")
    print("Derived Time-Averaged Stored Energy in the Magnetic Field (W_H):")
    sp.pretty_print(W_H)
    
    # The code below is for comparing with the options, not part of the main derivation output.
    # Let's construct the expressions from option D to compare.
    common_denominator = 2 * (omega/c) * (n**2 - 1) * ((n**2 + 1)*sin_theta**2 - 1) * sp.sqrt(n**2*sin_theta**2 - 1)
    
    W_E_D = (n**2 * (2*n**2*sin_theta**2 - 1)) / common_denominator * epsilon_0 * E_x0_i_amp**2
    W_H_D = (n**2 * (n**2*sin_theta**2 - 1)) / common_denominator * epsilon_0 * E_x0_i_amp**2
    
    print("\n" + "="*50 + "\n")
    print("Comparison with Option D:")
    if sp.simplify(W_E - W_E_D) == 0:
        print("The derived electric field energy W_E matches option D.")
    else:
        print("The derived electric field energy W_E does NOT match option D.")
        
    if sp.simplify(W_H - W_H_D) == 0:
        print("The derived magnetic field energy W_H matches option D.")
    else:
        print("The derived magnetic field energy W_H does NOT match option D.")
        print("The derived W_H is missing the factor (n^2*sin(theta)^2 - 1) in the numerator compared to option D.")

solve_evanescent_energy()