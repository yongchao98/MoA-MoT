import math
import numpy as np

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a roll of toilet paper to fall and completely unroll.

    The solution uses a 4th-order Runge-Kutta (RK4) method to solve the
    differential equations of motion for the roll, accounting for its changing
    mass, radius, and moment of inertia as it unrolls.
    """
    # 1. Define given constants
    G = 9.8  # Gravitational acceleration in m/s^2
    R_INNER = 0.04 / 2  # Radius of the cardboard cylinder in m
    T_PAPER = 0.5 / 1000  # Thickness of the paper in m
    N_WRAPS = 100  # Number of times the paper is wrapped
    M_P = 200 / 1000  # Mass of the paper in kg
    M_C = 20 / 1000  # Mass of the cardboard cylinder in kg

    # 2. Calculate derived constants
    # The initial outer radius of the full roll
    R_OUTER = R_INNER + N_WRAPS * T_PAPER

    # The total length of the paper, derived from the geometry.
    # The side-area of the paper annulus is pi*(R_OUTER^2 - R_INNER^2).
    # This area is also equal to L_TOTAL * T_PAPER.
    L_TOTAL = math.pi * (R_OUTER**2 - R_INNER**2) / T_PAPER

    # Moment of inertia of the cardboard tube, modeled as a thin shell (I = MR^2)
    I_C = M_C * R_INNER**2

    # Pre-calculate a common term for efficiency
    R0_SQR_MINUS_RI_SQR = R_OUTER**2 - R_INNER**2

    def get_acceleration(y):
        """
        Calculates the roll's instantaneous acceleration given the distance fallen y.
        """
        # Ensure y does not exceed the total length due to numerical precision
        if y >= L_TOTAL:
            y = L_TOTAL

        # Calculate the square of the current outer radius 'r' based on y
        # y = pi * (R_OUTER^2 - r^2) / T_PAPER  => r^2 = R_OUTER^2 - y*T_PAPER/pi
        r_sq = R_OUTER**2 - y * T_PAPER / math.pi
        
        # If the radius is less than the inner core, it means the paper has run out.
        if r_sq < R_INNER**2:
            r_sq = R_INNER**2

        # At the very end, the system is just the cardboard tube.
        # The acceleration simplifies to g / 2 (for a thin shell tube I=MR^2).
        if abs(r_sq - R_INNER**2) < 1e-9:
            return G * M_C / (M_C + I_C / R_INNER**2)

        # Calculate instantaneous mass of the remaining paper
        m_p_current = M_P * (r_sq - R_INNER**2) / R0_SQR_MINUS_RI_SQR
        m_current = M_C + m_p_current

        # Calculate instantaneous moment of inertia of the remaining paper
        # I_p = 1/2 * m_p * (r_outer^2 + r_inner^2)
        I_p_current = 0.5 * M_P * ((r_sq**2 - R_INNER**4) / R0_SQR_MINUS_RI_SQR)
        I_current = I_C + I_p_current

        # Calculate acceleration a = m*g / (m + I/r^2)
        acceleration = (m_current * G) / (m_current + I_current / r_sq)
        return acceleration

    def derivatives(t, state):
        """
        Calculates the derivatives [dy/dt, dv/dt] for the RK4 solver.
        """
        y, v = state
        acceleration = get_acceleration(y)
        return np.array([v, acceleration])

    # 3. Perform RK4 numerical simulation
    t = 0.0
    # State is a vector [y (position), v (velocity)]
    state = np.array([0.0, 0.0])
    dt = 0.001  # Time step in seconds for the simulation

    # Loop until the fallen distance 'y' exceeds the total paper length
    while state[0] < L_TOTAL:
        prev_state = state
        k1 = derivatives(t, state)
        k2 = derivatives(t + dt/2, state + dt/2 * k1)
        k3 = derivatives(t + dt/2, state + dt/2 * k2)
        k4 = derivatives(t + dt, state + dt * k3)
        state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + dt

    # 4. Interpolate to find the exact time of unrolling
    y_prev, v_prev = prev_state
    t_prev = t - dt

    # The remaining distance to fall is (L_TOTAL - y_prev)
    # Time to cover this distance is approximately distance/velocity
    time_to_go = (L_TOTAL - y_prev) / v_prev
    t_final = t_prev + time_to_go

    print("The toilet paper roll becomes fully unrolled after falling a distance of {:.2f} meters.".format(L_TOTAL))
    print("\nTo find the precise time, we look at the last simulation step:")
    print("Equation: time = t_previous + (L_total - y_previous) / v_previous")
    print("Numbers:  time = {:.4f} s + ({:.4f} m - {:.4f} m) / {:.4f} m/s".format(t_prev, L_TOTAL, y_prev, v_prev))
    print("\nThis gives the final time, rounded to the nearest hundredth of a second:")
    print("Time to unroll: {:.2f} seconds".format(t_final))

    # Final answer in the required format
    print("\n<<<" + "{:.2f}".format(t_final) + ">>>")

solve_toilet_paper_fall()