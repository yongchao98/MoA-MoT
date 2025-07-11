import math
from scipy.optimize import brentq

def solve_moment_capacity():
    """
    This function determines the moment at collapse for the given non-standard section.
    It follows the principles of limit state design for reinforced concrete.
    """
    # Step 1: Define Material Properties and Design Strengths
    # Given material characteristic strengths
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    # Partial safety factors (EC2)
    gamma_c = 1.5
    gamma_s = 1.15
    # Coefficients for concrete stress block
    alpha_cc = 0.85
    lambda_param = 0.8
    # Ultimate strain for concrete and steel properties
    ecu = 0.0035
    E_s = 200000.0  # MPa

    # Calculate design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s

    # Step 2: Define Section Geometry and Reinforcement
    # Geometry (all units in mm)
    b_top = 100.0
    b_bottom = 400.0
    h_trap = 300.0
    h_total = 400.0

    # Reinforcement bars
    d_bar = 20.0
    A_s_bar = math.pi * (d_bar / 2)**2
    
    # Steel layers are defined by their depth from the top compression fiber
    steel_layers = [
        {'d': 50.0, 'area': 2 * A_s_bar, 'id': 1},   # Top layer
        {'d': 260.0, 'area': 2 * A_s_bar, 'id': 2},  # Middle layer
        {'d': 350.0, 'area': 3 * A_s_bar, 'id': 3},  # Bottom layer
    ]

    # Helper function to find the section width at any depth y from the top
    def get_width(y):
        if y <= h_trap:
            return b_top + (b_bottom - b_top) * (y / h_trap)
        else:
            return b_bottom

    # Step 3: Find Neutral Axis Depth (x) by Enforcing Force Equilibrium
    # The equilibrium function calculates the net axial force for a given neutral axis depth 'x'
    def equilibrium(x):
        s = lambda_param * x  # Depth of the rectangular stress block

        # Calculate the compressive force in the concrete (Fc)
        if s <= h_trap: # Stress block is within the trapezoidal part
            Ac = (b_top + get_width(s)) / 2.0 * s
        else: # Stress block extends into the rectangular part
            Ac_trap = (b_top + b_bottom) / 2.0 * h_trap
            Ac_rect = b_bottom * (s - h_trap)
            Ac = Ac_trap + Ac_rect
        Fc = fcd * Ac
        
        net_force = Fc
        
        # Calculate forces in each steel layer and add to the net force
        for layer in steel_layers:
            d_i, A_si = layer['d'], layer['area']
            if d_i < x:  # Steel is in compression
                strain = ecu * (x - d_i) / x
                stress = min(fyd, E_s * strain) # Stress cannot exceed yield strength
                force = A_si * stress
                # Subtract the force of displaced concrete if the bar is within the stress block
                if d_i <= s:
                    force -= A_si * fcd
                net_force += force
            else:  # Steel is in tension
                strain = ecu * (d_i - x) / x
                stress = min(fyd, E_s * strain)
                force = A_si * stress
                net_force -= force
        return net_force

    # Use a numerical solver to find the root of the equilibrium function, which is the value of x
    try:
        x_sol = brentq(equilibrium, 0.1, h_total)
    except ValueError:
        print("Error: Could not find a solution for the neutral axis depth.")
        return 0

    # Step 4: Calculate the Moment at Collapse (M_rd)
    # With the neutral axis 'x' known, we can calculate all forces and their lever arms.
    s = lambda_param * x_sol

    # Calculate concrete force (Fc) and the location of its centroid (zc)
    if s <= h_trap:
        width_at_s = get_width(s)
        Ac = (b_top + width_at_s) / 2.0 * s
        zc = (s / 3.0) * (b_top + 2 * width_at_s) / (b_top + width_at_s)
    else:
        Ac_trap = (b_top + b_bottom) / 2.0 * h_trap
        zc_trap = (h_trap / 3.0) * (b_top + 2 * b_bottom) / (b_top + b_bottom)
        Ac_rect = b_bottom * (s - h_trap)
        zc_rect = h_trap + (s - h_trap) / 2.0
        Ac = Ac_trap + Ac_rect
        zc = (Ac_trap * zc_trap + Ac_rect * zc_rect) / Ac
    Fc = fcd * Ac

    # Calculate forces in each steel layer
    forces = {}
    for layer in steel_layers:
        d_i, A_si, id = layer['d'], layer['area'], layer['id']
        if d_i < x_sol:  # Compression
            strain = ecu * (x_sol - d_i) / x_sol
            stress = min(fyd, E_s * strain)
            if d_i <= s:
                force = A_si * (stress - fcd)
            else:
                force = A_si * stress
            forces[id] = {'type': 'C', 'val': force, 'd': d_i}
        else:  # Tension
            strain = ecu * (d_i - x_sol) / x_sol
            stress = min(fyd, E_s * strain)
            force = A_si * stress
            forces[id] = {'type': 'T', 'val': force, 'd': d_i}
    
    # Calculate total compression (C) and tension (T) forces and their centroids (z_C, z_T)
    C_total = Fc
    C_moment_about_top = Fc * zc
    T_total = 0
    T_moment_about_top = 0

    for id, f in forces.items():
        if f['type'] == 'C':
            C_total += f['val']
            C_moment_about_top += f['val'] * f['d']
        else:
            T_total += f['val']
            T_moment_about_top += f['val'] * f['d']
            
    z_C = C_moment_about_top / C_total if C_total != 0 else 0
    z_T = T_moment_about_top / T_total if T_total != 0 else 0

    # The lever arm is the distance between the resultant tension and compression forces
    lever_arm = z_T - z_C

    # The moment at collapse is the total tension force multiplied by the lever arm
    moment_rd_Nmm = T_total * lever_arm
    moment_rd_kNm = moment_rd_Nmm / 1e6

    # Step 5: Output the Final Equation and Result
    print("The moment at collapse (M_rd) is calculated as the total tension force (T) multiplied by the internal lever arm (z).")
    print(f"M_rd = T * z")
    print(f"M_rd = {T_total / 1000:.2f} kN * {lever_arm / 1000:.3f} m")
    print(f"M_rd = {moment_rd_kNm:.2f} kNm")
    
    return moment_rd_kNm

if __name__ == '__main__':
    final_answer = solve_moment_capacity()
    print(f"<<<{final_answer:.2f}>>>")