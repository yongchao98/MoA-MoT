import math

def solve_puzzle():
    """
    This function solves the multi-part physics and geometry puzzle to find the
    number of triangular prisms that fit in a cube.
    """
    
    # Part 1: Calculate T from heat transfer problem
    # ------------------------------------------------
    
    # Given parameters for the solar collector
    L = 1.5      # m, Length
    B = 0.85     # m, Width
    U_inf = 1.0  # m/s, Wind speed
    
    # Given air properties at film temperature
    nu_f = 15.11e-6 # m^2/s, Kinematic viscosity
    k_f = 0.0257    # W/(m*K), Thermal conductivity
    Pr_f = 0.707    # Prandtl number

    print("### Part 1: Calculation of Prism Base T ###")

    # The surface temperature theta_w(x) varies along the length.
    # To find the total heat loss, we use the average surface temperature.
    # Average of theta_w(x) = 30 + 10*sin(pi*x/L) over x=[0,L] is 30 + 10*(2/pi).
    theta_w_avg = 30 + 20 / math.pi
    print(f"1. Average surface temperature: theta_w_avg = 30 + 20/pi = {theta_w_avg:.3f} C")

    # The ambient temperature theta_infinity(y) varies with height.
    # Average of theta_inf(y) = 10 + 0.05y over y=[0,B] is 10 + 0.025*B.
    theta_inf_avg = 10 + 0.025 * B
    print(f"2. Average ambient temperature: theta_inf_avg = 10 + 0.025*{B} = {theta_inf_avg:.3f} C")

    # Average temperature difference drives the heat transfer
    delta_T_avg = theta_w_avg - theta_inf_avg
    print(f"3. Average temperature difference: Delta_T = {theta_w_avg:.3f} - {theta_inf_avg:.3f} = {delta_T_avg:.3f} K")

    # Calculate Reynolds number to determine the flow regime (laminar or turbulent)
    Re_L = (U_inf * L) / nu_f
    print(f"4. Reynolds number: Re_L = ({U_inf} * {L}) / {nu_f:.2e} = {Re_L:.0f}")
    # Since Re_L < 5e5, the flow is laminar.

    # Calculate the average Nusselt number for laminar flow over a flat plate
    Nu_L_avg = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    print(f"5. Average Nusselt number: Nu_L_avg = 0.664 * {Re_L:.0f}**0.5 * {Pr_f}**(1/3) = {Nu_L_avg:.3f}")

    # Calculate the average heat transfer coefficient
    h_avg = (Nu_L_avg * k_f) / L
    print(f"6. Average heat transfer coefficient: h_avg = ({Nu_L_avg:.3f} * {k_f}) / {L} = {h_avg:.3f} W/(m^2*K)")
    
    # Calculate total convective heat loss
    Area = L * B
    Q_dot_V = h_avg * Area * delta_T_avg
    print(f"7. Total heat loss: Q_dot_V = {h_avg:.3f} * {Area:.3f} * {delta_T_avg:.3f} = {Q_dot_V:.3f} W")

    # Calculate T and round to the nearest integer
    T_float = Q_dot_V / 80.0
    T = int(round(T_float))
    print(f"8. Prism base T = round({Q_dot_V:.3f} / 80) = round({T_float:.3f}) = {T}\n")

    # Part 2: Calculate D from beam mechanics problem
    # -----------------------------------------------
    print("### Part 2: Calculation of Prism Depth D ###")

    # Given parameters for the beam
    q0 = 3.0 # N/m, uniformly distributed load
    l = 2.0  # m, beam length
    
    # The maximum normal stress is sigma = M_max * z_max / I_y.
    # For a simply supported beam with uniform load, M_max = q0*l^2/8.
    # From the problem's definition of 'a', we deduce I_y = a, and z_max = 2a.
    # This simplifies sigma_max to (q0*l^2/8 * 2a) / a = q0*l^2/4.
    sigma_xx_max = (q0 * l**2) / 4.0
    print(f"1. Maximum normal stress: sigma_xx_max = ({q0} * {l}**2) / 4 = {sigma_xx_max:.3f} N/m^2")

    # Calculate D
    D_float = sigma_xx_max / 3.0
    D = int(D_float)
    print(f"2. Prism depth D = {sigma_xx_max:.3f} / 3 = {D_float:.3f} = {D}\n")
    
    # Part 3: Calculate how many prisms fit in the cube
    # ----------------------------------------------------
    print("### Part 3: Packing Prisms in the Cube ###")
    
    cube_dim = 3
    
    # A right triangular prism with legs T x T and depth D fits in a T x T x D bounding box.
    # Two such prisms fit perfectly into one bounding box.
    print(f"The prism's bounding box is {T}x{T}x{D}. The cube is {cube_dim}x{cube_dim}x{cube_dim}.")

    # Number of boxes that fit along each dimension of the cube
    num_x = math.floor(cube_dim / T)
    num_y = math.floor(cube_dim / T)
    num_z = math.floor(cube_dim / D)

    # Total number of boxes that can be packed
    num_boxes = num_x * num_y * num_z
    
    # Each box contains 2 prisms
    total_prisms = 2 * num_boxes
    
    print("\nFinal Calculation:")
    print(f"Number of prisms = 2 * (Number of boxes)")
    print(f"                 = 2 * (floor({cube_dim}/{T}) * floor({cube_dim}/{T}) * floor({cube_dim}/{D}))")
    print(f"                 = 2 * ({num_x} * {num_y} * {num_z})")
    print(f"                 = 2 * {num_boxes}")
    print(f"                 = {total_prisms}")

solve_puzzle()
<<<54>>>