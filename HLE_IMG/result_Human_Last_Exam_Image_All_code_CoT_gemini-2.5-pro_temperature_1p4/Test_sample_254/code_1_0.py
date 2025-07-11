import math

def solve_satellite_power():
    """
    Calculates the power incident on a photovoltaic cell on the Moon's surface
    after light from a source is reflected by two satellites.
    """
    # Constants and given values in SI units
    P = 1.0 * 10**9      # Power of the source in Watts (1 GW)
    S = 10.0             # Area of the photovoltaic cell in m^2
    M = 7.35 * 10**22    # Mass of the Moon in kg
    R = 1738 * 10**3     # Radius of the Moon in meters (1738 km)
    G = 6.67 * 10**-11   # Gravitational constant in m^3 kg^-1 s^-2
    T = 12 * 3600        # Orbital period in seconds (12 hours)
    pi = math.pi

    # Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law.
    # For the symmetric configuration required by the problem, the orbit must be circular,
    # so the orbital radius 'r' is equal to the semi-major axis 'a'.
    # T^2 / a^3 = 4 * pi^2 / (G * M) => a = (G * M * T^2 / (4 * pi^2))^(1/3)
    r = (G * M * T**2 / (4 * pi**2))**(1/3)

    # Step 2: Determine the effective distance using the virtual source method.
    # The geometric analysis shows that at the moment of reflection, the points O (Moon's center),
    # A, X, Y, B are collinear. The reflections at mirrors X and Y are at normal incidence.
    # We can track the light as if it comes from a series of virtual sources.
    # Let the origin be O. A is at R, X is at r, Y is at -r, B is at -R.
    # 1st virtual source A' (image of A in mirror X): A' = 2*r - R
    # 2nd virtual source A'' (image of A' in mirror Y): A'' = 2*(-r) - (2*r - R) = -4*r + R
    # The distance from the final virtual source A'' to the target B is:
    # d_A''B = |B - A''| = |(-R) - (-4*r + R)| = 4*r - 2*R
    d_effective = 4 * r - 2 * R

    # Step 3: Calculate the final power P' incident on the cell.
    # The original power P radiates into a hemisphere (2*pi sr). The final power on the cell S
    # is given by the formula P' = (P * S) / (2 * pi * d_effective^2).
    P_prime = (P * S) / (2 * pi * d_effective**2)

    # Step 4: Convert to microwatts for the final answer.
    P_prime_microwatts = P_prime * 1e6

    # --- Output the calculation steps ---
    print("--- Calculation Steps ---")
    print(f"Given source power P = {P:.1e} W")
    print(f"Given cell area S = {S:.1f} m^2")
    print(f"Given lunar mass M = {M:.2e} kg")
    print(f"Given lunar radius R = {R:.3e} m")
    print(f"Given orbital period T = {T} s")

    print("\n1. Calculate orbital radius 'r' (from Kepler's 3rd Law):")
    print(f"r = (G * M * T^2 / (4 * pi^2))^(1/3)")
    print(f"r = ({G:.2e} * {M:.2e} * {T}^2 / (4 * {pi:.4f}^2))^(1/3) = {r:.4e} m")
    
    print("\n2. Determine the effective distance from the final virtual source to the target:")
    print("The geometric configuration is collinear: A-X-O-Y-B.")
    print(f"d_effective = 4 * r - 2 * R")
    print(f"d_effective = 4 * {r:.4e} m - 2 * {R:.4e} m = {d_effective:.4e} m")
    
    print("\n3. Calculate the power P' incident on the cell:")
    print(f"P' = (P * S) / (2 * pi * d_effective^2)")
    print(f"P' = ({P:.1e} W * {S:.1f} m^2) / (2 * {pi:.4f} * ({d_effective:.4e} m)^2)")
    print(f"P' = {P_prime:.4e} W")
    
    print("\n4. Convert the result to microwatts:")
    print(f"P' (microwatts) = {P_prime:.4e} W * 1,000,000 = {P_prime_microwatts:.1f} uW")
    
    # Return the final answer in the specified format
    return f"<<<{P_prime_microwatts:.1f}>>>"

# Run the solver and print the final answer
final_answer = solve_satellite_power()
print(f"\n{final_answer}")
