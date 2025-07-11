import sympy as sp

def calculate_pyramid_flux():
    """
    Calculates the flux of the vector field F through the two yellow sides of a square pyramid.
    """
    # Define symbolic variables
    u, v = sp.symbols('u v')
    x, y, z = sp.symbols('x y z')

    # Define the vector field F
    F = sp.Matrix([3*x**3*y**2*z, 3*x**2*y**3, z])
    
    # --- Side 1: Yellow side in the positive y direction ---
    # Vertices: (-1, 1, 0), (1, 1, 0), (0, 0, 4)
    # Parameterization r1(u, v) for the surface
    x1 = (1 - v) * (2 * u - 1)
    y1 = 1 - v
    z1 = 4 * v
    r1 = sp.Matrix([x1, y1, z1])

    # Calculate the normal vector
    r1_u = r1.diff(u)
    r1_v = r1.diff(v)
    N1_raw = r1_u.cross(r1_v)
    
    # The normal vector should point outwards. For this face, the y-component should be positive.
    # N1_raw has a negative y-component, so we flip its sign.
    N1 = -N1_raw

    # Substitute the parameterization into the vector field F
    F_on_S1 = F.subs([(x, x1), (y, y1), (z, z1)])

    # Calculate the dot product F . N
    integrand1 = F_on_S1.dot(N1)
    
    # Integrate over the parameter domain u=[0,1], v=[0,1]
    flow1 = sp.integrate(integrand1, (u, 0, 1), (v, 0, 1))

    # --- Side 2: Yellow side in the negative y direction ---
    # Vertices: (-1, -1, 0), (1, -1, 0), (0, 0, 4)
    # Parameterization r2(u, v) for the surface
    x2 = (1 - v) * (2 * u - 1)
    y2 = -(1 - v)
    z2 = 4 * v
    r2 = sp.Matrix([x2, y2, z2])

    # Calculate the normal vector
    r2_u = r2.diff(u)
    r2_v = r2.diff(v)
    N2 = r2_u.cross(r2_v)

    # The normal vector should point outwards. For this face, the y-component should be negative.
    # N2 has a negative y-component, which is correct for this side.
    
    # Substitute the parameterization into the vector field F
    F_on_S2 = F.subs([(x, x2), (y, y2), (z, z2)])

    # Calculate the dot product F . N
    integrand2 = F_on_S2.dot(N2)
    
    # Integrate over the parameter domain u=[0,1], v=[0,1]
    flow2 = sp.integrate(integrand2, (u, 0, 1), (v, 0, 1))
    
    # --- Calculate and Print Total Flow ---
    total_flow = flow1 + flow2
    
    print(f"The energy flow through the first yellow side (y > 0) is: {flow1}")
    print(f"The energy flow through the second yellow side (y < 0) is: {flow2}")
    print("\nThe total energy flow through both yellow sides is the sum of the flows through each side.")
    print(f"Total Flow = {flow1} + {flow2} = {total_flow}")

if __name__ == '__main__':
    calculate_pyramid_flux()
    # The final answer as a simplified fraction is 104/21
    final_answer_float = 104/21
    print(f"\nThe final numerical answer is {final_answer_float:.4f}")
    print(f"<<<{104/21}>>>")
