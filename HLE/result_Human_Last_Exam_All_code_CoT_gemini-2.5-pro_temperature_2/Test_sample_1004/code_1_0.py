import sympy

def solve_boat_problem():
    """
    Calculates the final downstream distance of the boat based on the derived formula.
    
    The problem can be solved by integrating the river's flow velocity with respect
    to time for the two parts of the boat's journey: outward and return.
    
    1. Define river flow v_flow(y) as a piecewise function.
    2. The time to cross is determined by the boat's speed `v` (dt = dy/v).
    3. The downstream drift `x` is the integral of v_flow(y) * dt.
    
    - Drift outward (y from 0 to 3L/4): x1 = Integral(v_flow(y)/v dy) = (7*v0*L)/(16*v)
    - Drift return (y from 3L/4 to 0): x2 = Integral(v_flow(y)/v dy) = (7*v0*L)/(16*v)
    - Total drift = x1 + x2 = (7*v0*L)/(8*v)

    This code implements the final derived formula.
    """

    # --- Define variables (using example values) ---
    # L: width of the river (meters)
    # v0: maximum flow velocity of the river at the center (m/s)
    # v: speed of the boat relative to the water (m/s)
    
    # We use symbolic variables first to display the formula
    L_sym, v0_sym, v_sym = sympy.symbols('L v_0 v')

    # Example numerical values for calculation
    L = 160  # meters
    v0 = 4   # m/s
    v = 2    # m/s

    # The derived formula has constants 7 and 8
    numerator_const = 7
    denominator_const = 8

    # --- Calculation ---
    distance = (numerator_const * v0 * L) / (denominator_const * v)

    # --- Output Results ---
    print("This script calculates the final distance between the boat's start and end points.")
    print("-" * 70)
    print("The final derived formula for the distance is:")
    print(f"Distance = ({numerator_const} * {v0_sym} * {L_sym}) / ({denominator_const} * {v_sym})\n")
    
    print("Using the following example values:")
    print(f"L (River Width) = {L} m")
    print(f"v_0 (Max Flow Velocity) = {v0} m/s")
    print(f"v (Boat Speed) = {v} m/s\n")
    
    print("Plugging the values into the formula:")
    print(f"Distance = ({numerator_const} * {v0} * {L}) / ({denominator_const} * {v})")
    
    print("\nFinal calculated distance:")
    print(f"{distance} meters")


# Run the solver
solve_boat_problem()