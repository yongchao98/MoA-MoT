import math

def solve_period():
    """
    This function calculates and prints the formula for the period of the pendulum-like system.
    
    The derivation steps are as follows:
    1.  The kinetic energy (T) and potential energy (U) of the system are determined.
        -   Coordinates: x (horizontal position of top disk's center), θ (rod's angle with horizontal).
        -   Mass of each disk: m. Radius: R. Rod length: 4R. Gravity: g.
        -   Kinetic Energy of Disk 1 (rolling): T1 = (5/8)*m*v_x^2
        -   Kinetic Energy of Disk 2: T2 = T_trans + T_rot
            -   T_trans = (1/2)*m*((v_x - 4*R*dθ*sin(θ))^2 + (4*R*dθ*cos(θ))^2)
            -   T_rot = (1/8)*m*v_x^2 (since ω2 = ω1 = v_x/R)
        -   Total T = (5/4)*m*v_x^2 - 4*m*R*v_x*dθ*sin(θ) + 8*m*R^2*dθ^2
        -   Potential Energy U = mgR + mg(R + 4*R*sin(θ)) = 2*mgR + 4*mgR*sin(θ)
    2.  The Lagrangian L = T - U is linearized for small oscillations around the equilibrium position (θ = -π/2).
        Let θ = -π/2 + φ, where φ is a small angle.
    3.  The linearized equations of motion are derived from the Lagrangian.
        (5/2)*m*d^2x + 4*m*R*d^2φ = 0
        4*m*R*d^2x + 16*m*R^2*d^2φ + 4*m*g*R*φ = 0
    4.  Solving these coupled equations gives the angular frequency ω for the oscillatory mode:
        ω^2 = (5*g) / (12*R)
    5.  The period T_period is calculated as 2*π / ω.
        T_period = 2*π * sqrt(12*R / (5*g))
    """
    
    # The numbers in the final equation for the period
    num_two = 2
    num_twelve = 12
    num_five = 5
    
    # Print the final formula in a descriptive way
    print("The derived period of motion, T, is given by the formula:")
    print(f"T = {num_two} * pi * sqrt(({num_twelve} * R) / ({num_five} * g))")

solve_period()