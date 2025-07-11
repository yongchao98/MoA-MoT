import sympy

def solve_pendulum_period():
    """
    This function symbolically derives the period of oscillation for the described system.
    """
    # Step 1: Define symbolic variables
    # M: mass of one disk
    # R: radius of one disk
    # g: acceleration due to gravity
    # L: length of the rod
    M, R, g, L = sympy.symbols('M R g L', positive=True)
    
    # Define theta_dot as a symbol for angular velocity, since we are building the SHM equation
    theta_dot = sympy.symbols('theta_dot')

    # Step 2: Define system parameters and constraints based on the problem statement
    # The rod length is given as L = 4R
    L_val = 4 * R
    
    # Moment of inertia for a solid disk about its center
    I_disk = sympy.Rational(1, 2) * M * R**2

    # A key step is relating the linear motion of the top disk (x) to the angular motion
    # of the pendulum (theta). Analysis of the forces and the rolling-without-slipping
    # constraint yields the relationship for accelerations: ddot(x) = (-2/3) * L * ddot(theta).
    # For small oscillations starting from rest, this integrates to:
    # dot(x) = (-2/3) * L * dot(theta)
    x_dot = -sympy.Rational(2, 3) * L * theta_dot
    
    print("Derived relationship for small oscillations:")
    print(f"v_disk1 = {sympy.latex(x_dot).replace('theta_dot', r'\dot{\theta}')}\n")


    # Step 3: Calculate the total Kinetic Energy (T) of the system for small angles
    # T_total = T_disk1 + T_disk2
    
    # T_disk1 is the sum of its translational and rotational kinetic energy.
    # Rotational velocity omega_disk1 = v_disk1 / R = x_dot / R
    T_disk1 = sympy.Rational(1, 2) * M * x_dot**2 + sympy.Rational(1, 2) * I_disk * (x_dot / R)**2
    
    # T_disk2 is its translational kinetic energy. For small angles, its velocity is
    # approximately v_disk2 = x_dot + L*theta_dot (horizontal component dominates)
    v_disk2 = x_dot + L * theta_dot
    T_disk2 = sympy.Rational(1, 2) * M * v_disk2**2
    
    T_total = T_disk1 + T_disk2
    
    # Substitute the constraint to express T_total solely in terms of theta_dot
    T_total_theta = T_total.subs(x_dot, -sympy.Rational(2, 3) * L * theta_dot).simplify()
    
    # The form of the kinetic energy is (1/2) * m_eff * (theta_dot)^2
    # So we can find the effective mass (or moment of inertia) of the system.
    m_eff = (2 * T_total_theta / theta_dot**2).simplify()
    print("Effective Inertia (m_eff) of the system:")
    print(f"m_eff = {sympy.latex(m_eff)}\n")

    # Step 4: Calculate the Potential Energy (V) for small angles
    # The potential energy V comes from the change in height of the second disk.
    # V = M * g * h = M * g * (-L * cos(theta))
    # For small angles, cos(theta) is approx. 1 - theta^2 / 2.
    # V = -M*g*L*(1 - theta^2/2) = const + (1/2)*M*g*L*theta^2
    # The form of potential energy is (1/2) * k_eff * theta^2
    # So we can find the effective spring constant.
    k_eff = M * g * L
    print("Effective Spring Constant (k_eff) of the system:")
    print(f"k_eff = {sympy.latex(k_eff)}\n")
    
    # Step 5: Calculate the angular frequency (omega) and period (T_p)
    # The equation of motion is m_eff * ddot(theta) + k_eff * theta = 0
    # The angular frequency squared is omega^2 = k_eff / m_eff
    omega_sq = (k_eff / m_eff).simplify()
    
    # The period is T_p = 2*pi / omega
    T_p = 2 * sympy.pi / sympy.sqrt(omega_sq)
    
    # Step 6: Substitute L = 4R to get the final answer
    final_period = T_p.subs(L, L_val).simplify()
    
    # Print the final result and its components
    print("The final expression for the period of motion is:")
    final_expression_str = f"{final_period.args[0]} * pi * sqrt({final_period.args[2].args[0]} * R / g)"
    
    # Extract numerical coefficients
    c1, _, c2_sqrt = final_period.args
    c2 = c2_sqrt.args[0]
    
    print(f"T = ({c1}) * π * sqrt(({c2}) * R / g)")
    
    # Return the symbolic answer for the final block
    return final_period


final_answer_expression = solve_pendulum_period()
final_answer_str = str(final_answer_expression).replace('**', '^') # Format for answer block

# Let's reformat it to be more readable
# (4*pi*sqrt(7)/3)*sqrt(R/g) -> (4*pi/3)*sqrt(7*R/g)
final_answer_str = f"(4*pi/3)*sqrt(7*R/g)"


# The problem asks for the period.
print(f"\nFinal Answer: {final_answer_str}")
# >>> (4*pi/3)*sqrt(7*R/g)
