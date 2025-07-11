import math
import numpy as np

def solve_toilet_paper_problem():
    """
    Solves the falling toilet paper problem using a Runge-Kutta (RK4) numerical method.
    """
    # 1. Define physical constants of the system
    g = 9.81  # Acceleration due to gravity (m/s^2)
    r_c = 0.02  # Radius of the cardboard cylinder (m)
    m_c = 0.02  # Mass of the cardboard cylinder (kg)
    M_p_initial = 0.2  # Initial mass of the paper (kg)
    paper_thickness = 0.0005  # Thickness of the paper (m)
    N_wraps = 100  # Number of wraps

    # Calculate derived constants
    # Initial outer radius of the full roll
    R_initial = r_c + N_wraps * paper_thickness
    # Total length of the paper using the volume conservation method
    L_total = (math.pi * (R_initial**2 - r_c**2)) / paper_thickness
    # Moment of inertia of the cylinder (assuming a thin hollow tube)
    I_c = m_c * r_c**2

    # 2. Define system properties as functions of unrolled length 'y'
    def R_sq(y):
        """Returns the square of the outer radius as a function of unrolled length y."""
        # Clamp y to avoid math domain errors from overshooting
        y_clamped = min(y, L_total)
        # The area unrolled (y*t) equals the change in cross-sectional area
        val = R_initial**2 - y_clamped * paper_thickness / math.pi
        return max(val, r_c**2) # Radius cannot be smaller than the cylinder's radius

    def R(y):
        """Returns the outer radius."""
        return math.sqrt(R_sq(y))

    def M_p(y):
        """Returns the remaining mass of the paper."""
        y_clamped = min(y, L_total)
        return M_p_initial * (1 - y_clamped / L_total)

    def M(y):
        """Returns the total mass of the falling roll (cylinder + paper)."""
        return m_c + M_p(y)

    def I(y):
        """Returns the moment of inertia of the falling roll."""
        # I_paper for a hollow cylinder is 0.5 * M * (R_outer^2 + R_inner^2)
        i_paper = 0.5 * M_p(y) * (R_sq(y) + r_c**2)
        return I_c + i_paper

    def dR_dy(y):
        """Returns the derivative of the radius with respect to y."""
        current_R = R(y)
        # If paper is gone, radius is constant.
        if current_R <= r_c + 1e-9:
            return 0.0
        # From R^2 = C - k*y, we get 2R*dR/dy = -k => dR/dy = -k/(2R)
        return -paper_thickness / (2 * math.pi * current_R)

    def dI_dy(y):
        """Returns the derivative of the moment of inertia with respect to y."""
        if y >= L_total:
            return 0.0
        # Use product rule: d/dy [ M_p(y) * (R_sq(y) + r_c^2) ]
        dMp_dy = -M_p_initial / L_total
        dR2_dy = -paper_thickness / math.pi
        term1 = dMp_dy * (R_sq(y) + r_c**2)
        term2 = M_p(y) * dR2_dy
        return 0.5 * (term1 + term2)

    def acceleration(y, v):
        """
        Calculates the roll's acceleration a = f(y, v) using the full dynamic equations.
        a = ( M*g + v^2*[dM/dy - d(I/R)/dy] ) / ( M + I/R^2 ) -> Expanded form below.
        """
        if y >= L_total:
            return 0.0
        
        current_M = M(y)
        current_R = R(y)
        current_I = I(y)
        current_dR_dy = dR_dy(y)
        current_dI_dy = dI_dy(y)
        
        # This term comes from d(Mv)/dt = Ma + v*dM/dt
        # We use dM/dy = -M_p_initial / L_total
        mass_change_term_force = -(M_p_initial / L_total) * v**2

        # The term with dI/dt and dR/dt in the torque equation gives this corrective term for tension
        tension_correction_force = - ( (v**2 * current_dI_dy / current_R) + 
                                     (current_I * v * (-v / current_R**2 * current_dR_dy)) )

        # Denominator of tension T = I*a/R^2 + correction
        # Net force Mg-T = Ma_eff -> Mg - I*a/R^2 - correction = Ma_eff
        
        numerator = current_M * g + mass_change_term_force - tension_correction_force
        denominator = current_M + current_I / current_R**2
        
        if denominator == 0:
            return 0.0
        return numerator / denominator
        
    # Simplified acceleration function from expanded formula.
    def full_acceleration(y,v):
        if y >= L_total: return 0.0
        
        cM = M(y); cR = R(y); cI = I(y); cdRdy = dR_dy(y); cdIdy = dI_dy(y)

        if cR < r_c + 1e-9: return 0.0
        
        v_sq_term_bracket = ((M_p_initial / L_total) 
                             - (cdIdy / cR**2) 
                             + (cI * cdRdy / cR**3))
        
        numerator = cM * g + v**2 * v_sq_term_bracket
        denominator = cM + cI / cR**2
        if denominator == 0: return 0
        return numerator / denominator
        
    def dS_dt(state):
        y, v = state
        a = full_acceleration(y, v)
        return np.array([v, a])

    # 3. RK4 solver implementation
    t = 0.0
    dt = 0.002  # Time step in seconds for accuracy
    state = np.array([0.0, 0.0])  # Initial state [y, v]
    
    # Store history to interpolate at the end
    y_history = [0.0]
    t_history = [0.0]

    while state[0] < L_total:
        k1 = dS_dt(state)
        k2 = dS_dt(state + 0.5 * dt * k1)
        k3 = dS_dt(state + 0.5 * dt * k2)
        k4 = dS_dt(state + dt * k3)
        
        state += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t += dt
        
        t_history.append(t)
        y_history.append(state[0])

    # 4. Interpolate for the final exact time
    y_last = y_history[-2]
    t_last = t_history[-2]
    y_final = y_history[-1]
    t_final_approx = t_history[-1]

    # Linear interpolation: t = t1 + (t2-t1) * (y_target-y1) / (y2-y1)
    t_actual = t_last + (t_final_approx - t_last) * (L_total - y_last) / (y_final - y_last)
    
    # Print the requested information and final answer
    print("This problem is solved by numerically integrating the equations of motion.")
    print("The final time is found by interpolating the last simulation step to find the exact time when the unrolled length equals the total paper length.")
    print("\nConstants:")
    print(f"Total paper length (L_total): {L_total:.4f} m")

    print("\nFinal interpolation step:")
    print(f"Time (t) when unrolled length ({y_last:.4f} m) was just under the total: {t_last:.4f} s")
    print(f"Time (t) when unrolled length ({y_final:.4f} m) was just over the total: {t_final_approx:.4f} s")
    print(f"\nThe equation for the final interpolated time is:")
    print(f"t_final = {t_last:.4f} + ({t_final_approx:.4f} - {t_last:.4f}) * ({L_total:.4f} - {y_last:.4f}) / ({y_final:.4f} - {y_last:.4f})")
    
    print(f"\nCalculated final time to unroll completely:")
    print(f"{t_actual:.2f} seconds")

    return t_actual

# Run the simulation and capture the final answer
final_time = solve_toilet_paper_problem()
print(f"<<<{final_time:.2f}>>>")