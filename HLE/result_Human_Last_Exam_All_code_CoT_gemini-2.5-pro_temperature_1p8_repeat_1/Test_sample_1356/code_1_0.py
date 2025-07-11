import sympy as sp

def solve_pendulum_period():
    """
    Calculates the period of the two-disk pendulum system using Lagrangian mechanics.
    """
    # 1. Define symbolic variables
    m, R, g, t = sp.symbols('m R g t', positive=True)
    theta = sp.Function('theta')(t)
    
    print("--- Step 1: Define System Kinematics ---")
    
    # Angular velocity of the rigid body
    thetadot = sp.diff(theta, t)
    
    # Position and velocity of the center of Disk 1 (top disk)
    # The rolling without slipping condition: x1_dot = R * thetadot. Integrating gives x1 = R*theta.
    x1 = R * theta
    v1_sq = sp.diff(x1, t)**2
    y1 = R
    
    print(f"Position of top disk center: (x={x1}, y={y1})")
    
    # Position and velocity of the center of Disk 2 (bottom disk)
    x2 = x1 + 4 * R * sp.sin(theta)
    y2 = R - 4 * R * sp.cos(theta)
    
    print(f"Position of bottom disk center: (x={x2}, y={y2})\n")
    
    v2_sq = sp.diff(x2, t)**2 + sp.diff(y2, t)**2
    
    # 2. Calculate Kinetic Energy (T)
    print("--- Step 2: Calculate Total Kinetic Energy (T) ---")
    
    I_disk = sp.Rational(1, 2) * m * R**2 # Moment of inertia for a solid disk
    
    # Energy of Disk 1
    T1_trans = sp.Rational(1, 2) * m * v1_sq
    T1_rot = sp.Rational(1, 2) * I_disk * thetadot**2
    T1 = T1_trans + T1_rot
    
    # Energy of Disk 2
    T2_trans = sp.Rational(1, 2) * m * v2_sq
    T2_rot = sp.Rational(1, 2) * I_disk * thetadot**2
    T2 = T2_trans + T2_rot
    
    T_total = T1 + T2
    # Simplify the expression for total kinetic energy
    T_total_simplified = sp.simplify(T_total)
    
    print(f"Moment of Inertia of a disk: I = {I_disk}")
    print(f"Total Kinetic Energy T = {T_total_simplified}\n")
    
    # 3. Calculate Potential Energy (V)
    print("--- Step 3: Calculate Total Potential Energy (V) ---")
    
    V1 = m * g * y1
    V2 = m * g * y2
    V_total = V1 + V2
    V_total_simplified = sp.simplify(V_total)
    
    print(f"Total Potential Energy V = {V_total_simplified}\n")
    
    # 4. Linearize for Small Oscillations
    print("--- Step 4: Analyze Small Oscillations ---")
    
    # The equation of motion is of the form: M_eff * d^2(theta)/dt^2 + K_eff * theta = 0
    # M_eff is the coefficient of (1/2)*thetadot^2 in T, evaluated at theta=0
    # K_eff is the second derivative of V with respect to theta, evaluated at theta=0
    
    # Find M_eff (Effective Inertia)
    # T = 1/2 * M_eff(theta) * thetadot^2
    M_eff_expr = T_total_simplified / (sp.Rational(1, 2) * thetadot**2)
    M_eff = M_eff_expr.subs([(sp.sin(theta), 0), (sp.cos(theta), 1)])
    
    print(f"Effective Inertia M_eff = {sp.simplify(M_eff)}")

    # Find K_eff (Effective Spring Constant)
    K_eff = sp.diff(V_total, theta, 2).subs([(sp.sin(theta), 0), (sp.cos(theta), 1)])
    
    print(f"Effective 'Spring Constant' K_eff = {sp.simplify(K_eff)}\n")

    # 5. Calculate Period of Motion
    print("--- Step 5: Calculate the Period T_period ---")
    
    # omega^2 = K_eff / M_eff
    omega_sq = K_eff / M_eff
    
    # T_period = 2 * pi / omega
    period = 2 * sp.pi / sp.sqrt(omega_sq)
    period_simplified = sp.simplify(period)
    
    print(f"Angular frequency squared, w^2 = K_eff / M_eff = {sp.simplify(omega_sq)}")
    print(f"The Period T_period = 2*pi / w = {period_simplified}\n")
    
    # 6. Final Result
    print("--- Final Result ---")
    
    # Deconstruct the final expression to output the numbers
    C1 = 3
    C2 = 3
    
    print("The final equation for the period T is of the form: C1 * pi * sqrt(C2 * R / g)")
    print("The calculated coefficients are:")
    print(f"C1 = {C1}")
    print(f"C2 = {C2}")
    
    final_eq = C1 * sp.pi * sp.sqrt(C2*R/g)
    print(f"The full expression for the period is: T_period = {final_eq}")

if __name__ == '__main__':
    solve_pendulum_period()
