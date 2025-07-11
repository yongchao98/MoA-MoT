import math

def calculate_force():
    """
    Calculates the instantaneous force on a coil in a cylindrical magnetic structure.
    """
    # Define physical parameters with example values
    R = 0.1  # Radius of the outer coil in meters
    N = 100  # Turns of the second (AC) coil
    N0 = 500  # Turns of the first (DC) coil
    I0 = 10.0  # DC current in Amperes
    i0 = 0.5  # Amplitude of AC current in Amperes
    omega = 100 * math.pi  # Angular frequency in rad/s
    g = 0.005  # Radial gap in meters (5 mm)
    mu0 = 4 * math.pi * 1e-7  # Permeability of free space in H/m
    alpha_T = 0.002  # Temperature coefficient of permeability in 1/K
    T = 350.0  # Operating temperature in Kelvin (T1)
    T0 = 300.0  # Reference temperature in Kelvin
    Bs = 1.5  # Saturation flux density in Tesla
    t = 0.0025  # Time in seconds

    # Step 1: Calculate time-varying current i(t)
    i_t = i0 * math.sin(omega * t)

    # Step 2: Calculate the numerator of the force equation fraction
    # This combines mu0, temperature correction, coil parameters, and currents.
    temp_factor = (1 - alpha_T * (T - T0))
    numerator = mu0 * temp_factor * N0 * I0 * i_t

    # Step 3: Calculate the denominator of the force equation fraction
    # This includes the gap and the saturation effect.
    saturation_term = 1 + (mu0 * N0 * I0) / (g * Bs)
    denominator = g**2 * saturation_term

    # Step 4: Calculate the final instantaneous force
    force = -2 * math.pi * R * N * (numerator / denominator)
    
    # --- Output the results ---
    print("Derived Force Equation (from Choice B):")
    print("f_x(t) = -2*pi*R*N * [mu0*(1-alpha_T*(T-T0))*N0*I0*i(t)] / [g^2*(1 + (mu0*N0*I0)/(g*Bs))]\n")
    
    print(f"Calculating the force at t = {t} s with the following values:")
    print(f"   R (radius) = {R} m")
    print(f"   N (AC turns) = {N}")
    print(f"   N0 (DC turns) = {N0}")
    print(f"   I0 (DC current) = {I0} A")
    print(f"   i0 (AC current amplitude) = {i0} A")
    print(f"   omega (frequency) = {omega:.2f} rad/s")
    print(f"   g (gap) = {g} m")
    print(f"   mu0 (permeability) = {mu0:.3e} H/m")
    print(f"   alpha_T (temp. coeff) = {alpha_T} 1/K")
    print(f"   T (temperature) = {T} K")
    print(f"   T0 (ref temperature) = {T0} K")
    print(f"   Bs (saturation flux) = {Bs} T")
    print("-" * 30)

    # Print the calculated values for each part of the equation
    print("Calculated components of the equation:")
    print(f"   i(t) = {i0}*sin({omega:.2f}*{t}) = {i_t:.4f} A")
    print(f"   Numerator (inside fraction): {numerator:.3e}")
    print(f"   Denominator (inside fraction): {denominator:.3e}")
    print(f"   Constant Multiplier (-2*pi*R*N): {-2*math.pi*R*N:.2f}")
    print("-" * 30)
    print(f"Final instantaneous force f_x({t}) = {force:.4f} N")

if __name__ == "__main__":
    calculate_force()
