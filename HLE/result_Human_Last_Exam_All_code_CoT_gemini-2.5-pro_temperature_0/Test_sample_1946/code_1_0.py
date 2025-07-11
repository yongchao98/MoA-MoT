import numpy as np
from scipy.integrate import solve_ivp

def solve_toilet_paper_fall():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.
    """
    # 1. Define physical constants from the problem in SI units
    g = 9.81  # Acceleration due to gravity (m/s^2)
    d_cylinder_cm = 4
    paper_thickness_mm = 0.5
    n_wraps = 100
    m_paper_g = 200
    m_cylinder_g = 20
    h_initial_m = 50 # This is a distractor, as the paper length determines the fall distance

    # Convert to SI units
    r_c = (d_cylinder_cm / 100) / 2  # m (radius of cardboard cylinder)
    d_paper = paper_thickness_mm / 1000  # m (thickness of one paper layer)
    m_p = m_paper_g / 1000  # kg (mass of paper)
    m_c = m_cylinder_g / 1000  # kg (mass of cardboard cylinder)

    # 2. Derive secondary constants
    r_0 = r_c + n_wraps * d_paper  # Initial outer radius of the roll
    
    # Total length of the paper, L = pi * (r_0^2 - r_c^2) / d_paper
    l_total = np.pi * (r_0**2 - r_c**2) / d_paper

    # Pre-calculate a constant term for efficiency in the model
    r0_sq_minus_rc_sq = r_0**2 - r_c**2

    # 3. Define the system of ODEs: model(t, state)
    # The state vector is [y, v, r]
    def model(t, state):
        y, v, r = state
        
        # Protect against the radius becoming smaller than the cardboard core
        if r <= r_c:
            # This case is handled by the event function, but as a safeguard:
            # The roll is empty, it would fall at acceleration g.
            dv_dt = g
            dr_dt = 0.0
        else:
            # Calculate instantaneous mass of the roll
            m_paper_inst = m_p * (r**2 - r_c**2) / r0_sq_minus_rc_sq
            m_roll_inst = m_c + m_paper_inst
            
            # Calculate instantaneous moment of inertia
            # I_cylinder (thin wall approx) = m_c * r_c^2
            # I_paper (hollow cylinder) = 0.5 * m_paper_inst * (r^2 + r_c^2)
            I_inst = m_c * r_c**2 + 0.5 * m_paper_inst * (r**2 + r_c**2)
            
            # Calculate linear acceleration using a = g / (1 + I / (m * r^2))
            dv_dt = g / (1 + I_inst / (m_roll_inst * r**2))
            
            # Calculate the rate of change of the radius
            # v = - (2*pi*r / d_paper) * dr/dt  =>  dr/dt = -v * d_paper / (2*pi*r)
            dr_dt = -v * d_paper / (2 * np.pi * r) if v > 0 else 0.0

        dy_dt = v
        
        return [dy_dt, dv_dt, dr_dt]

    # 4. Define the event to stop the simulation when the roll is empty
    def end_of_roll(t, state):
        # This function will be zero when y equals the total length
        return state[0] - l_total
    end_of_roll.terminal = True  # Stop the integration when this event occurs
    end_of_roll.direction = 1    # Trigger when event function goes from - to +

    # 5. Run the simulation
    initial_state = [0, 0, r_0]  # y=0, v=0, r=r_0
    t_span = [0, 15]  # A reasonable time span for the simulation

    solution = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=initial_state,
        events=end_of_roll,
        dense_output=True
    )

    # 6. Output the results
    print("--- Problem Parameters ---")
    print(f"Cardboard inner diameter: {d_cylinder_cm} cm")
    print(f"Paper thickness: {paper_thickness_mm} mm")
    print(f"Number of wraps: {n_wraps}")
    print(f"Total paper mass: {m_paper_g} g")
    print(f"Cardboard mass: {m_cylinder_g} g")
    
    print("\n--- Derived Parameters (SI Units) ---")
    print(f"Gravity (g): {g} m/s^2")
    print(f"Cardboard radius (r_c): {r_c:.4f} m")
    print(f"Paper thickness (d): {d_paper:.4f} m")
    print(f"Paper mass (M_p): {m_p:.3f} kg")
    print(f"Cardboard mass (M_c): {m_c:.3f} kg")
    print(f"Initial roll radius (r_0): {r_0:.4f} m")
    print(f"Total paper length (L_total): {l_total:.2f} m")
    print("------------------------------------")

    if solution.status == 1:  # Status 1 means an event was triggered
        time_to_unroll = solution.t_events[0][0]
        print(f"\nRESULT:")
        print(f"The toilet paper takes {time_to_unroll:.2f} seconds to unroll completely.")
        
        # As requested, printing the numbers in the "final equation"
        # can be interpreted as the final state of the system.
        final_state = solution.sol(time_to_unroll)
        y_final, v_final, r_final = final_state
        print("\n--- Final State ---")
        print(f"Time (t): {time_to_unroll:.2f} s")
        print(f"Distance fallen (y): {y_final:.2f} m")
        print(f"Velocity (v): {v_final:.2f} m/s")
        print(f"Roll radius (r): {r_final:.4f} m (final) vs {r_c:.4f} m (core)")
        
        return time_to_unroll
    else:
        print("\nSolver did not find the event. Check parameters or increase t_span.")
        print(f"Solver status: {solution.status}, Message: {solution.message}")
        return None

if __name__ == '__main__':
    final_time = solve_toilet_paper_fall()
    if final_time is not None:
        # The final answer is printed here in the required format.
        # The code above provides the detailed breakdown.
        pass # The final answer will be appended outside the code block.