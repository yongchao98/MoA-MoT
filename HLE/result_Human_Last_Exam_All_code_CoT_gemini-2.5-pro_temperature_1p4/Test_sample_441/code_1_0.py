import numpy as np

def solve_hipparchus_problem():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    approximating the angular motion of an object moving on a square.
    """

    # Part 1: Determine phi (frequency ratio)
    # The motion has 4-fold symmetry. The deviation from uniform angular motion
    # completes 4 cycles in one full orbit.
    # The deferent frequency is omega_d (1 cycle/orbit).
    # The deviation frequency in the model is (omega_e - omega_d).
    # Equating them: omega_e - omega_d = 4 * omega_d
    # This gives omega_e = 5 * omega_d.
    phi = 5.0

    # Part 2: Determine R (radii ratio)
    # R is found by matching the amplitude of the angular deviation.
    # The model's max deviation is A_model = arcsin(1/R).
    # We must find the actual max deviation A_actual.

    # Let the square have vertices at (+-1, +-1). The perimeter is 8.
    # Let the object move at speed v=1, so the period T=8.
    # The mean angular velocity (deferent frequency) is omega_d = 2*pi/T = pi/4.

    # The extrema of the angular deviation occur when the instantaneous
    # angular velocity equals the mean angular velocity.
    # On the side from (1,-1) to (1,1), parameterized by arclength s in [0,2],
    # the angular velocity d(theta)/ds = 1 / (1 + (s-1)^2).
    # With v=1, ds/dt=1, so d(theta)/dt = d(theta)/ds.
    # Set d(theta)/dt = omega_d:
    # 1 / (1 + (s-1)^2) = pi/4
    # (s-1)^2 = 4/pi - 1
    # s = 1 +- sqrt(4/pi - 1)
    
    # Let u = sqrt(4/pi - 1). The deviation amplitude is given by the formula:
    # A_actual = arctan(u) - (pi/4)*u
    
    val_4_div_pi = 4 / np.pi
    u = np.sqrt(val_4_div_pi - 1)
    
    term1 = np.arctan(u)
    term2 = (np.pi / 4) * u
    
    # The amplitude of the angular deviation in radians
    A_actual = term1 - term2
    
    # Now, solve for R: A_actual = arcsin(1/R) => R = 1 / sin(A_actual)
    sin_A = np.sin(A_actual)
    R = 1 / sin_A

    # Output the results, including intermediate steps for R's calculation.
    print("The ordered pair is (R, phi).")
    print(f"phi = {phi}")
    print("\nCalculation for R:")
    print("R = 1 / sin(A_actual)")
    print("A_actual = arctan(u) - (pi/4)*u")
    print("u = sqrt(4/pi - 1)")
    print("-" * 20)
    print(f"u = sqrt({val_4_div_pi:.6f} - 1) = sqrt({u**2:.6f}) = {u:.6f}")
    print(f"A_actual = arctan({u:.6f}) - (pi/4)*{u:.6f}")
    print(f"A_actual = {term1:.6f} - {term2:.6f} = {A_actual:.6f} radians")
    print(f"R = 1 / sin({A_actual:.6f}) = 1 / {sin_A:.6f}")
    print(f"R = {R:.4f}")
    print("-" * 20)
    print(f"The calculated ordered pair is ({R:.4f}, {phi}).")


solve_hipparchus_problem()
<<<({14.0338}, {5.0})>>>