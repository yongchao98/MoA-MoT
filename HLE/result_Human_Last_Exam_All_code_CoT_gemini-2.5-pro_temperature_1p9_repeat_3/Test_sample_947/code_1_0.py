import sympy as sp

def solve_superconductor_field():
    """
    Calculates and prints the symbolic expression for the magnetic field
    around a stack of superconducting strips for |x| >> a.
    """

    # --- Step 1: Define all symbolic variables ---
    # These represent the physical quantities in the problem.
    x, z = sp.symbols('x z', real=True)              # Position coordinates
    Ha = sp.Symbol('H_a', positive=True)             # Applied magnetic field
    Jc = sp.Symbol('J_c', positive=True)             # Critical current density
    d = sp.Symbol('d', positive=True)                # Strip thickness
    w = sp.Symbol('w', positive=True)                # Strip half-width
    D = sp.Symbol('D', positive=True)                # Stacking interval
    pi = sp.pi

    # --- Step 2: Define intermediate and given parameters ---
    # Define the characteristic field H0 and the penetration depth 'a'.
    # H0 is given in the problem statement.
    H0 = Jc * d / pi

    # The relationship between the applied field Ha and the flux penetration
    # depth 'a' for a single thin strip is given by Ha = H0 * arccosh(w/a).
    # We solve this for 'a'. This approximation ignores the fields from
    # neighboring strips when determining 'a'.
    a = w / sp.cosh(Ha / H0)

    # --- Step 3: Model the induced currents and calculate the field ---
    # In the far-field limit (|x| >> w), each strip's induced current
    # can be modeled as two opposite line currents.
    # The total current on one side of a strip is I = (width of current) * (sheet current density)
    I = (w - a) * (Jc * d)

    # The field from an infinite stack of such line-current pairs can be
    # calculated using the complex potential method. We define a constant C
    # that encapsulates the source strength.
    C = (I / (2 * D)) * sp.sinh(pi * (w + a) / D)

    # Define dimensionless coordinates u and v for simplicity.
    u = pi * x / D
    v = pi * z / D

    # --- Step 4: Construct the expressions for the field components ---
    # The field H is a vector with H_x and H_z components. These components
    # are derived from the complex potential of the current array.
    
    # Common denominator for both H_x and H_z components.
    denominator = (sp.sinh(u)**2 + sp.sin(v)**2)**2

    # Numerator for the induced H_x component.
    H_x_ind_numerator = C * sp.Rational(1, 2) * (sp.cosh(2*u) * sp.cos(2*v) - 1)

    # Numerator for the induced H_z component.
    H_z_ind_numerator = C * sp.Rational(1, 2) * sp.sinh(2*u) * sp.sin(2*v)

    # The total H_x field is just the induced field.
    H_x = H_x_ind_numerator / denominator
    
    # The total H_z field is the sum of the applied field and the induced field.
    H_z_induced = H_z_ind_numerator / denominator
    H_z = Ha + H_z_induced

    # --- Step 5: Print the final expressions ---
    # We display the symbolic expressions for the components of the magnetic field.
    # To make it more readable, we print the definitions of a, C, u, and v first.

    print("The solution is based on a far-field approximation (|x| >> w), where each strip's current is modeled as a magnetic dipole.")
    print("The penetration depth 'a' is determined by the applied field Ha:")
    sp.pprint(sp.Eq(sp.Symbol('a'), a), use_unicode=True)
    print("\nThe field components are expressed in terms of the following constants and variables:")
    sp.pprint(sp.Eq(sp.Symbol('C'), C, evaluate=False), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('u'), u), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('v'), v), use_unicode=True)
    
    print("\nThe x-component of the magnetic field H_x(x, z) is:")
    sp.pprint(sp.Eq(sp.Function('H_x')(x, z), H_x), use_unicode=True)
    
    print("\nThe z-component of the magnetic field H_z(x, z) is:")
    sp.pprint(sp.Eq(sp.Function('H_z')(x, z), H_z), use_unicode=True)

solve_superconductor_field()