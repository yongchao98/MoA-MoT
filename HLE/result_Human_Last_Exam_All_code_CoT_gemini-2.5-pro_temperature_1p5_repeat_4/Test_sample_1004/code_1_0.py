import sympy

def solve_boat_problem():
    """
    Calculates and explains the solution to the river boat problem.
    """
    # Define symbolic variables
    L, v0, v, y = sympy.symbols('L v_0 v y')

    # 1. Model the river's flow velocity
    # The velocity profile is parabolic: v_flow(y) = k*y*(L-y)
    # At y=L/2, v_flow(L/2) = v0.
    # v0 = k * (L/2) * (L - L/2) = k * L**2 / 4  =>  k = 4*v0 / L**2
    v_flow = (4 * v0 / L**2) * y * (L - y)

    # 2. Calculate downstream drift for the outbound trip (from y=0 to y=3L/4)
    # The time to cross dy is dt = dy / v
    # The drift dx is v_flow(y) * dt
    # The total drift is the integral of (v_flow/v) dy
    integrand_out = v_flow / v
    # The boat travels to a distance L/4 from the opposite bank, so y goes from 0 to L - L/4 = 3L/4
    delta_x_out = sympy.integrate(integrand_out, (y, 0, 3 * L / 4))

    # 3. Calculate downstream drift for the return trip (from y=3L/4 to y=0)
    # The physics is the same, so the drift over the same y-interval is identical.
    delta_x_in = delta_x_out

    # 4. Calculate the total displacement
    total_distance = delta_x_out + delta_x_in

    # 5. Print the results
    print("Step-by-step derivation of the boat's downstream displacement:")
    print("-" * 60)
    
    print("1. The river's flow velocity profile is modeled as:")
    print(f"   v_flow(y) = {v_flow}\n")

    print("2. The downstream displacement during the outbound trip (from y=0 to y=3L/4) is:")
    print(f"   Δx_out = Integral({integrand_out}) dy from 0 to 3L/4")
    print(f"   Δx_out = {delta_x_out}\n")
    
    print("3. The downstream displacement during the return trip (from y=3L/4 to y=0) is the same:")
    print(f"   Δx_in = {delta_x_in}\n")

    print("4. The total distance is the sum of the two displacements:")
    print(f"   Total Distance = Δx_out + Δx_in = {total_distance}\n")
    
    # Extracting the numerical coefficients for the final output as requested
    numer, denom = total_distance.as_numer_denom()
    numer_coeff = numer / (v0 * L)
    denom_coeff = denom / v

    print("-" * 60)
    print("Final Answer:")
    print("The final equation for the total downstream distance is:")
    print(f"    ({int(numer_coeff)} * v_0 * L)")
    print(f"    ---------")
    print(f"    ({int(denom_coeff)} * v)")

if __name__ == '__main__':
    solve_boat_problem()