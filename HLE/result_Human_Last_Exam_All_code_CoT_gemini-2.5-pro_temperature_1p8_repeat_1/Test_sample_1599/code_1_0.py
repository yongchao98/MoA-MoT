import numpy as np
from scipy.optimize import fsolve

def solve_acceleration():
    """
    Calculates the man's final acceleration based on the problem statement.
    """
    # Given constants
    u = 1.5  # m/s
    t0 = 0
    t1 = 4
    a1 = -0.15  # m/s^2
    t3 = 15
    a2 = 0.25  # m/s^2
    t5 = 23
    t6 = 40
    alpha_deg = 130
    gamma_deg = 40

    # Step 1: Solve for intermediate times t2 and t4
    # t2 is when the man's speed becomes 0 after decelerating.
    # u_t2 = u_t1 + a1 * (t2 - t1) => 0 = 1.5 - 0.15 * (t2 - 4)
    t2 = u / abs(a1) + t1

    # t4 is when the man's speed reaches u again after accelerating from rest.
    # u_t4 = u_t3 + a2 * (t4 - t3) => 1.5 = 0 + 0.25 * (t4 - 15)
    t4 = u / a2 + t3

    # Step 2: Calculate man's position at key time points
    y_m_t1 = u * (t1 - t0)
    y_m_t2 = y_m_t1 + u * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    y_m_t3 = y_m_t2
    y_m_t4 = y_m_t3 + 0.5 * a2 * (t4 - t3)**2  # Starts from rest at t3
    y_m_t5 = y_m_t4 + u * (t5 - t4)

    # Convert angles to radians
    alpha = np.deg2rad(alpha_deg)
    gamma = np.deg2rad(gamma_deg)

    # Step 3: Solve for bird's speed 'v' and initial ground speed 'vg'
    # This system is derived from the "planned rendezvous" conditions.
    # It relies on the identity that sin(130)=cos(40) and cos(130)=-sin(40)
    def equations(p):
        v, vg = p
        if v**2 < vg**2 or v <= 0 or vg < 0:
            return [1e6, 1e6] # Invalid physical solution
        
        vz = np.sqrt(v**2 - vg**2)
        Cg, Sg, Tg = np.cos(gamma), np.sin(gamma), np.tan(gamma)
        
        # Bird's coordinates at t4 based on its segmented journey
        x_b_t4 = 4 * vg * Cg + 10 * v
        y_b_t4 = v - 4 * vg * Sg
        z_b_t4 = 4 * vz + 6 * v

        # Equation 1: From the planned dive geometry, z_b(t4) = x_b(t4) * tan(gamma)
        eq1 = z_b_t4 - (x_b_t4 * Tg)
        
        # Equation 2: From equating the bird's and man's positions at the planned rendezvous time
        # y_b(t4) must equal y_m(planned_t_f)
        planned_dt = x_b_t4 / (v * Cg) # This is (t_f_planned - t4)
        eq2 = y_b_t4 - (y_m_t4 + u * planned_dt)
        
        return [eq1, eq2]

    # Numerically solve the system of equations
    initial_guess = [30, 2] # A reasonable guess based on preliminary analysis
    v, vg = fsolve(equations, initial_guess)

    # Step 4: Calculate the bird's position at t5
    vz = np.sqrt(v**2 - vg**2)
    Cg, Sg = np.cos(gamma), np.sin(gamma)
    
    x_b_t4 = 4 * vg * Cg + 10 * v
    y_b_t4 = v - 4 * vg * Sg
    z_b_t4 = 4 * vz + 6 * v
    
    vx_plan = -v * Cg
    vz_plan = -v * Sg
    
    x_b_t5 = x_b_t4 + vx_plan * (t5 - t4)
    y_b_t5 = y_b_t4
    z_b_t5 = z_b_t4 + vz_plan * (t5 - t4)

    # Step 5: Calculate the final rendezvous position y_m(t6)
    dt_final = t6 - t5
    # The bird travels a distance of v*dt_final from P_b(t5) to P_final=(0, y_m(t6), 0)
    # (v*dt_final)^2 = (x_b_t5)^2 + (y_m_t6 - y_b_t5)^2 + (z_b_t5)^2
    # The new trajectory is "northward", so y_m(t6) must be greater than y_b(t5)
    dy_final_sq = (v * dt_final)**2 - x_b_t5**2 - z_b_t5**2
    y_m_t6 = y_b_t5 + np.sqrt(dy_final_sq)
    
    # Step 6: Use man's kinematics to solve for a3
    # y_m(t6) = y_m(t5) + u * dt_final + 0.5 * a3 * dt_final^2
    a3 = (y_m_t6 - y_m_t5 - u * dt_final) / (0.5 * dt_final**2)

    # --- Final Output ---
    print(f"The initial speed of the man is u = {u:.3f} m/s.")
    print(f"The man stops at t2 = {t2:.3f} s and starts accelerating again at t3 = {t3:.3f} s.")
    print(f"The man reaches speed u again at t4 = {t4:.3f} s.")
    print(f"The man's position at t5 = {t5:.3f} s is y = {y_m_t5:.3f} m.")
    print(f"By solving the rendezvous equations, the bird's constant speed is found to be v = {v:.3f} m/s.")
    print(f"At the moment of the wind gust (t5 = {t5:.3f} s), the bird is at y = {y_b_t5:.3f} m.")
    print(f"In the final leg, the bird and man meet at y = {y_m_t6:.3f} m at t6 = {t6:.3f} s.")
    print("\nTo find the man's final acceleration, we use the equation: y(t6) = y(t5) + u*(t6-t5) + 0.5*a3*(t6-t5)^2")
    final_dt = t6-t5
    final_dt_sq = final_dt**2
    y_term = y_m_t6
    y0_term = y_m_t5
    ut_term = u * final_dt
    print(f"Substituting the values:")
    print(f"{y_term:.3f} = {y0_term:.3f} + {u:.3f}*{final_dt:.3f} + 0.5*a3*{final_dt_sq:.3f}")
    print(f"{y_term:.3f} = {y0_term:.3f} + {ut_term:.3f} + {0.5*final_dt_sq:.3f}*a3")
    numerator = y_term - y0_term - ut_term
    denominator = 0.5 * final_dt_sq
    print(f"a3 = ({y_term:.3f} - {y0_term:.3f} - {ut_term:.3f}) / {denominator:.3f}")
    print(f"a3 = {numerator:.3f} / {denominator:.3f}")
    print(f"a3 = {a3:.3f} m/s^2")

solve_acceleration()