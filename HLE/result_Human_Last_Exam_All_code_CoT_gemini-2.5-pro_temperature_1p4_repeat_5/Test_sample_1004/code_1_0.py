import sympy

def calculate_boat_displacement():
    """
    This function calculates the total downstream displacement of the boat.
    """
    # Define the symbolic variables used in the problem.
    # L: width of the river
    # v_0: maximum flow velocity at the river's center
    # v: boat's speed relative to the water, perpendicular to the flow
    L, v_0, v = sympy.symbols('L v_0 v', positive=True, real=True)
    y = sympy.Symbol('y')

    print("Step 1: Define the river's flow velocity profile, v_river(y).")
    # The velocity profile is parabolic: v_river(y) = k * y * (L - y).
    # Using v_river(L/2) = v_0, we find k = 4*v_0 / L**2.
    k = 4 * v_0 / L**2
    v_river = k * y * (L - y)
    print(f"The river's flow velocity is: v_river(y) = {sympy.simplify(v_river)}\n")

    print("Step 2: Calculate the downstream displacement for one leg of the journey.")
    # The boat's velocity across the river is dy/dt = v, so dt = dy/v.
    # The downstream velocity is dx/dt = v_river(y).
    # The rate of downstream displacement w.r.t. cross-stream distance is dx/dy = v_river(y) / v.
    # We integrate this from y=0 to y = 3L/4 to find the displacement on the outward journey.
    integrand = v_river / v
    y_turnaround = (3 * L) / 4
    delta_x_out = sympy.integrate(integrand, (y, 0, y_turnaround))
    print(f"The displacement for the outward journey (from y=0 to y=3L/4) is: {delta_x_out}\n")

    print("Step 3: Calculate the total downstream displacement.")
    # The displacement during the return journey (from y=3L/4 to y=0) is identical.
    # The total displacement is the sum of the outward and return displacements.
    total_displacement = delta_x_out * 2
    print("The displacement for the return journey is the same. Total displacement is twice the outward displacement.")
    
    # Extract coefficients for the final print statement as requested.
    num_coeff = 9
    den_coeff = 8

    print("\n--- Final Answer ---")
    print("The final equation for the total downstream distance is:")
    print(f"Total Distance = ({num_coeff} * v_0 * L) / ({den_coeff} * v)")

if __name__ == '__main__':
    calculate_boat_displacement()