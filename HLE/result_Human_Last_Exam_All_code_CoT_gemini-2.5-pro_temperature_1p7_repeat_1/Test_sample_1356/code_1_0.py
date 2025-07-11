from fractions import Fraction

def solve_pendulum_period():
    """
    This function calculates the period of the two-disk pendulum system
    by deriving the linearized equation of motion using the Lagrangian method.
    """
    # Let M be the mass and R be the radius of each disk.
    # The Lagrangian L = T - V. We find the coefficients for the linearized E.O.M.
    # The form of the E.O.M is I_eff * theta_ddot + k_eff * theta = 0.
    # The period is then 2 * pi * sqrt(I_eff / k_eff).

    # 1. Kinetic Energy T = A*M*(dx/dt)^2 + B*M*R^2*(d(theta)/dt)^2 + C*M*R*(dx/dt)*(d(theta)/dt)
    # T_disk1 = 1/2*M*(dx/dt)^2 (trans) + 1/2*I*(omega_1)^2 (rot)
    #         = 1/2*M*(dx/dt)^2 + 1/2*(1/2*M*R^2)*(dx/dt / R)^2 = 3/4*M*(dx/dt)^2
    # T_disk2 = 1/2*M*v_2^2 + 1/2*I*(d(theta)/dt)^2
    # v_2^2 = (dx/dt + 4R*d(theta)/dt)^2 + (0)^2 for small angles
    # T_disk2 approx = 1/2*M*(dx/dt)^2 + 4*M*R*(dx/dt)*(d(theta)/dt) + 8*M*R^2*(d(theta)/dt)^2 + 1/4*M*R^2*(d(theta)/dt)^2
    
    # Coefficients for T/M, for small angles (cos(theta)=1)
    A = Fraction(3, 4) + Fraction(1, 2)            # for (dx/dt)^2
    B = Fraction(1, 2) * 16 + Fraction(1, 4)       # for R^2*(d(theta)/dt)^2
    C = Fraction(1, 2) * 8                         # for R*(dx/dt)*(d(theta)/dt)

    # 2. Potential Energy V = -4*M*g*R*cos(theta)
    # For small angles, cos(theta) approx 1 - theta^2/2.
    # V approx V_0 + 1/2 * k_eff * theta^2
    # V approx -4*M*g*R + 2*M*g*R*theta^2
    # k_eff is the coefficient of 1/2 * theta^2 in V, so k_eff = 4*M*g*R
    k_eff_coeff = Fraction(4)

    # 3. From dL/dx = 0, p_x = dL/d(x_dot) is conserved. System starts at rest, so p_x=0.
    # p_x = 2*A*M*(dx/dt) + C*M*R*(d(theta)/dt) = 0
    # => dx/dt = - (C / (2*A)) * R * d(theta)/dt
    dx_dt_relation_coeff = -C / (2 * A)

    # 4. Substitute into the theta equation of motion:
    # d/dt(dL/d(theta_dot)) - dL/d(theta) = 0
    # [2*B*M*R^2 - C^2/(2*A)*M*R^2] * theta_ddot + k_eff * theta = 0
    # I_eff = [2*B - C^2/(2*A)]*M*R^2
    
    I_eff_coeff = 2 * B - (C**2) / (2 * A)

    # 5. The angular frequency squared is omega^2 = k_eff / I_eff
    omega_sq_coeff_g_over_R = k_eff_coeff / I_eff_coeff

    # 6. The period P = 2*pi / omega = 2*pi * sqrt(1 / omega_sq_coeff_g_over_R * R/g)
    period_coeff_numerator = omega_sq_coeff_g_over_R.denominator
    period_coeff_denominator = omega_sq_coeff_g_over_R.numerator

    print("Step-by-step derivation of the period of motion:")
    print("-" * 50)
    print("1. The problem is solved using Lagrangian mechanics for small oscillations.")
    print("2. The equation of motion is simplified to the form: I_eff * d^2(theta)/dt^2 + k_eff * theta = 0.")
    print(f"3. The effective inertia term is found to be I_eff = ({I_eff_coeff}) * M * R^2.")
    print(f"4. The effective spring constant is found to be k_eff = ({k_eff_coeff}) * M * g * R.")
    print("5. The squared angular frequency is omega^2 = k_eff / I_eff.")
    print(f"   omega^2 = ({k_eff_coeff} * M*g*R) / ({I_eff_coeff} * M*R^2) = {omega_sq_coeff_g_over_R} * (g/R).")
    print("6. The period of motion is P = 2*pi / omega.")
    print("-" * 50)
    print("The final expression for the period is:")
    # The final print shows the numbers used in the equation as requested
    print(f"P = 2 * pi * sqrt( ({period_coeff_numerator} * R) / ({period_coeff_denominator} * g) )")

solve_pendulum_period()
<<<P = 2 * pi * sqrt( (101 * R) / (40 * g) )>>>