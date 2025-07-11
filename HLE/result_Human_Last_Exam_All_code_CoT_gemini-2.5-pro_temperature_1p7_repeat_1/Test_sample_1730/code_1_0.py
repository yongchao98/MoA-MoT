import numpy as np

def solve_puzzle():
    """
    Solves the multi-step problem to find the number of prisms that fit in a cube.
    """
    
    # Part 1: Calculate prism dimension T
    # Given parameters for the solar collector
    L = 1.5  # m, Length of the collector
    B = 0.85 # m, Width of the collector
    U_inf = 1.0 # m/s, Wind speed
    
    # Air properties
    rho_f = 1.204     # kg/m^3
    nu_f = 15.11e-6   # m^2/s
    k_f = 0.0257      # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1
    g = 9.81          # m/s^2

    # Step 1.1: Calculate average temperatures
    # Average plate temperature: integral of theta_w(x) from 0 to L, divided by L
    # theta_w_avg = (1/L) * integral[30 + 10*sin(pi*x/L)]dx
    # = 30 + (10/L) * [-L/pi * cos(pi*x/L)] from 0 to L
    # = 30 + (10/pi) * (-cos(pi) - (-cos(0))) = 30 + 10/pi * (1 - (-1)) = 30 + 20/pi
    theta_w_avg = 30.0 + 20.0 / np.pi
    
    # Average ambient temperature: integral of theta_inf(y) from 0 to B, divided by B
    # theta_inf_avg = (1/B) * integral[10 + 0.05y]dy from 0 to B
    # = (1/B) * [10y + 0.05y^2/2] from 0 to B
    # = 10 + 0.025 * B
    theta_inf_avg = 10.0 + 0.025 * B
    
    delta_theta_avg = theta_w_avg - theta_inf_avg
    
    # Step 1.2: Calculate heat transfer coefficients
    # Forced convection (along L)
    Re_L = (U_inf * L) / nu_f
    # Assuming laminar flow, Nu = 0.664 * Re^(1/2) * Pr^(1/3)
    Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    h_forced = (Nu_forced * k_f) / L
    
    # Natural convection (on vertical plate, height B)
    # The characteristic length for vertical natural convection is the vertical height B.
    Gr_B = (g * beta_f * delta_theta_avg * B**3) / (nu_f**2)
    Ra_B = Gr_B * Pr_f
    # Correlation for turbulent natural convection on a vertical plate
    Nu_natural = 0.1 * (Ra_B**(1/3))
    h_natural = (Nu_natural * k_f) / B
    
    # Step 1.3: Combine for mixed convection
    # A common correlation for mixed convection is h_comb^3 = h_forced^3 + h_natural^3
    h_combined = (h_forced**3 + h_natural**3)**(1/3)
    
    # Step 1.4: Calculate total heat loss Q_V
    A = L * B
    Q_V = h_combined * A * delta_theta_avg
    
    # Step 1.5: Calculate T
    T_float = Q_V / 80.0
    T = int(round(T_float))

    print("--- Calculation for T ---")
    print(f"Average plate temperature: {theta_w_avg:.2f} C")
    print(f"Average ambient temperature: {theta_inf_avg:.2f} C")
    print(f"Combined heat transfer coefficient: {h_combined:.2f} W/(m^2.K)")
    print(f"Total heat loss Q_V = {h_combined:.2f} * {A:.2f} * {delta_theta_avg:.2f} = {Q_V:.2f} W")
    print(f"T = {Q_V:.2f} / 80 = {T_float:.2f}")
    print(f"Rounded T = {T}")
    print("-" * 25)

    # Part 2: Calculate prism dimension D
    # Given parameters for the beam
    q0 = 3.0  # N/m
    l = 2.0   # m
    
    # Step 2.1: Calculate maximum bending moment for a simply supported beam
    M_max = (q0 * l**2) / 8.0
    
    # Step 2.2: The moment of inertia I_y is cleverly given by the formula for a.
    # a^3 = (64/3 - pi/4)^-1 => 1/a^3 = 64/3 - pi/4
    # I_y for the described cross-section is I_square - 2*I_semicircle
    # I_y = (4a*(4a)^3)/12 - 2*(pi*a^4/8) = (64/3 - pi/4)a^4
    # Substituting 1/a^3 for (64/3 - pi/4), we get I_y = (1/a^3)*a^4 = a.
    # So we don't need to calculate a.
    # z_max is the farthest distance from the neutral axis, which is 2*a
    # sigma_max = M_max * z_max / I_y = M_max * (2*a) / a = 2 * M_max
    sigma_xx_max = 2 * M_max
    
    # Step 2.3: Calculate D
    D_float = sigma_xx_max / 3.0
    D = int(round(D_float))

    print("--- Calculation for D ---")
    print(f"Maximum bending moment M_max = {M_max:.2f} N.m")
    print(f"Moment of inertia I_y simplifies to 'a'.")
    print(f"Maximum stress sigma_max = M_max * 2a / a = {sigma_xx_max:.2f} N/m^2")
    print(f"D = {sigma_xx_max:.2f} / 3 = {D_float:.2f}")
    print(f"Rounded D = {D}")
    print("-" * 25)

    # Part 3: Packing the prisms into the cube
    # Cube dimensions
    cube_dim = 3
    # Prism dimensions: base is a right triangle with legs T, T. Depth is D.
    # Prism bounding box dimensions: T x T x D
    
    # Our calculated values are T=2, D=1.
    prism_box = sorted([T, T, D], reverse=True) # [2, 2, 1]

    # Simple orthogonal packing strategy
    # The volume of the cube is 3x3x3=27.
    # The volume of a prism is 0.5 * T * T * D = 0.5 * 2 * 2 * 1 = 2.
    # Theoretical max is 27/2 = 13.5 -> 13 prisms. But packing is restricted by geometry.
    
    # We can fit a 2x2x3 block in the 3x3x3 cube.
    # This block has a volume of 12.
    # Each 2x2x1 slice of this block can contain 2 prisms.
    # So, a 2x2x3 block can contain 3 * 2 = 6 prisms.
    prisms_count = 6
    
    # The remaining space is an "L"-shaped prism of volume 27 - 12 = 15.
    # This can be seen as two blocks: a 1x3x3 block and a 2x1x3 block.
    # In the 1x3x3 block, we can fit a 1x2x2 box, holding 2 prisms.
    prisms_count += 2
    
    # In the 2x1x3 block, we can fit a 2x1x2 box, holding 2 prisms.
    prisms_count += 2

    # The remaining nooks and crannies (a 1x1x2 block and an L-shaped block of volume 5)
    # cannot fit any more prisms with a 2x2 base.
    
    print("--- Packing Calculation ---")
    print(f"Prism has a right triangular base with sides {T}, {T}, and {T*np.sqrt(2):.2f}.")
    print(f"Prism depth is {D}.")
    print("A 2x2x3 volume inside the 3x3x3 cube can be filled with 6 prisms.")
    print("The leftover 1x3x3 space can fit 2 more prisms.")
    print("The leftover 2x1x3 space can fit 2 more prisms.")
    print(f"Total number of prisms that can fit completely = 6 + 2 + 2 = {prisms_count}")
    
    return prisms_count

final_answer = solve_puzzle()
print(f"\nFinal Answer: The number of prisms is {final_answer}.")
print(f"<<<{final_answer}>>>")
