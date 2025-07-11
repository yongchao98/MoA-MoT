import sys
from io import StringIO

# A helper function to format the output
def print_derivation():
    """
    This function prints the step-by-step derivation of the solution.
    """
    print("### Derivation of the Downstream Distance ###\n")
    
    # Step 1: Model the river flow velocity
    print("Step 1: Define the River's Velocity Profile")
    print("Let the river width be L and the coordinate across the river be y (from 0 to L).")
    print("The river flow velocity v_flow(y) is:")
    print(" - Zero at the shores: v_flow(0) = 0 and v_flow(L) = 0.")
    print(" - Maximum (v0) at the center: v_flow(L/2) = v0.")
    print("A parabolic profile fits these conditions: v_flow(y) = k * y * (L - y).")
    print("Solving for k at the center: v0 = k * (L/2) * (L - L/2) => v0 = k * L^2/4 => k = 4*v0 / L^2.")
    print("So, the velocity profile is: v_flow(y) = (4*v0/L^2) * y * (L - y).\n")

    # Step 2: Calculate drift on the outbound trip
    print("Step 2: Calculate Downstream Drift on the Outbound Trip")
    print("The boat travels from y=0 to y = L - L/4 = 3L/4 with speed v perpendicular to the flow.")
    print("The boat's y-velocity is dy/dt = v, so we can write time as dt = dy/v.")
    print("The downstream drift (x_out) is the integral of the river's velocity over time:")
    print("x_out = ∫ v_flow(y) dt = ∫[0 to 3L/4] (4*v0/L^2)*y*(L-y) * (dy/v)")
    print("x_out = (4*v0)/(v*L^2) * ∫[0 to 3L/4] (L*y - y^2) dy")
    print("The integral of (L*y - y^2) is (L*y^2/2 - y^3/3).")
    print("Evaluating from 0 to 3L/4: (4*v0)/(v*L^2) * [L/2*(3L/4)^2 - 1/3*(3L/4)^3]")
    print("= (4*v0)/(v*L^2) * [9*L^3/32 - 9*L^3/64]")
    print("= (4*v0)/(v*L^2) * [9*L^3/64]")
    print("x_out = (36 * v0 * L) / (64 * v) = (9 * v0 * L) / (16 * v).\n")
    
    # Step 3: Calculate drift on the return trip
    print("Step 3: Calculate Downstream Drift on the Return Trip")
    print("The boat travels back from y=3L/4 to y=0.")
    print("The time taken and the river velocities experienced are the same as the outbound trip.")
    print("Therefore, the downstream drift on the return trip (x_in) is the same as x_out.")
    print("x_in = (9 * v0 * L) / (16 * v).\n")
    
    # Step 4: Calculate total distance
    print("Step 4: Calculate Total Downstream Distance")
    print("Total Distance D = x_out + x_in")
    print("D = (9*v0*L)/(16*v) + (9*v0*L)/(16*v)")
    print("D = 2 * (9*v0*L) / (16*v)")
    print("D = (18 * v0 * L) / (16 * v)\n")

    # Final Answer
    print("--- FINAL ANSWER ---")
    print("The total distance between the returning position and the original starting point is:")
    print("Distance = (9 * v0 * L) / (8 * v)")
    print("\nEach number in the final equation is:")
    print("Numerator: 9")
    print("Denominator: 8")

# Execute the function to print the solution
print_derivation()