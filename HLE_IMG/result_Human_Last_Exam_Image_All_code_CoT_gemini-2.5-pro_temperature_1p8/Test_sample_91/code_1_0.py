import math

def solve_drawbridge_velocity():
    """
    Calculates the vertical velocity of a drawbridge edge based on given parameters.
    """
    # Given parameters
    L = 40.0  # Length of the drawbridge in meters
    h = 10.0  # Height of the edge above the ground in meters

    # Step 1: Find cos(theta) and sin(theta) at the moment h = 10 m.
    # The geometric relationship is h = L * cos(theta).
    cos_theta = h / L
    # Using the identity sin^2(theta) + cos^2(theta) = 1
    sin_theta_val = math.sqrt(1 - cos_theta**2)

    # Step 2: Use the given expression for the rate of change of theta.
    # d(theta)/dt = -(3 * pi / 10) / cos(pi / 12)
    cos_pi_12 = math.cos(math.pi / 12)
    dtheta_dt_val = -(3 * math.pi / 10) / cos_pi_12

    # Step 3: Calculate the vertical velocity dh/dt.
    # The formula from differentiation is dh/dt = -L * sin(theta) * d(theta)/dt.
    dh_dt = -L * sin_theta_val * dtheta_dt_val

    # Print the explanation and the steps
    print("The vertical velocity of the drawbridge (dh/dt) is found from the height equation h = L * cos(θ).")
    print("Differentiating with respect to time, we get: dh/dt = -L * sin(θ) * dθ/dt\n")
    print("Step 1: Find sin(θ) when h = 10 m and L = 40 m.")
    print(f"cos(θ) = h / L = {h} / {L} = {cos_theta}")
    print(f"sin(θ) = sqrt(1 - cos²(θ)) = sqrt(1 - {cos_theta}²) = sqrt(15)/4\n")

    print("Step 2: State the given rate of change for θ.")
    print("dθ/dt = -(3 * π / 10) / cos(π / 12)\n")

    print("Step 3: Substitute these values into the velocity equation.")
    print("dh/dt = -L * sin(θ) * dθ/dt")
    print(f"dh/dt = -({L}) * (sqrt(15)/4) * [-(3 * π / 10) / cos(π/12)]")
    print("Simplifying the expression, we get:")
    print("dh/dt = (3 * π * sqrt(15)) / cos(π/12)\n")

    print("Final Calculation:")
    print(f"dh/dt = (3 * {math.pi:.4f} * {math.sqrt(15):.4f}) / {cos_pi_12:.4f}")
    print(f"The vertical velocity is approximately {dh_dt:.2f} m/s.")

solve_drawbridge_velocity()