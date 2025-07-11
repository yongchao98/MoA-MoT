import math

def solve_joule_heat():
    """
    This function explains and derives the formula for the Joule heat dissipated
    by a shrinking, charge-leaking sphere.
    """
    print("This problem requires calculating the total Joule heat dissipated by a sphere that is both shrinking and leaking charge.")
    print("The amount of heat depends on the exact process. We will use a physical model where leakage occurs at a constant critical electric field (E_crit) at the sphere's surface.\n")
    
    # Step 1: Define the model based on initial conditions
    print("--- Step 1: The Constant Electric Field Model ---")
    print("The electric field E at the surface of a sphere with radius r and potential V_sphere is E = V_sphere / r.")
    print("From the initial state (radius a, potential V), we define the critical field for leakage: E_crit = V / a.")
    print("We assume the process evolves such that the field at the changing surface r(t) is always E_crit.")
    print("    E_surface = V(r) / r = E_crit = V / a\n")

    # Step 2: Express charge and potential in terms of radius
    print("--- Step 2: Express Q(r) and V(r) ---")
    print("Using the model from Step 1, the potential at any radius r during the process is:")
    print("    V(r) = E_crit * r = (V/a) * r")
    print("The charge Q on a sphere is related to its potential by Q = 4 * pi * epsilon_0 * r * V(r). Substituting V(r):")
    print("    Q(r) = 4 * pi * epsilon_0 * r * (V/a) * r")
    print("    Q(r) = 4 * pi * epsilon_0 * (V/a) * r**2\n")

    # Step 3: Set up the integral for Joule Heat
    print("--- Step 3: Calculate the Joule Heat (H) ---")
    print("Joule heat is generated as charge leaks through the resistive atmosphere.")
    print("The incremental heat dH for a charge change dQ at potential V(r) is dH = -V(r) * dQ.")
    print("To find the total heat H, we integrate dH from the start (r=a) to the end (r=0).")
    print("First, we find dQ by differentiating Q(r) with respect to r:")
    print("    dQ = d/dr [4 * pi * epsilon_0 * (V/a) * r**2] dr")
    print("    dQ = 8 * pi * epsilon_0 * (V/a) * r * dr")
    print("Now, we set up the integral for H:")
    print("    H = integral from r=a to r=0 of [-V(r) * dQ]")
    print("    H = integral from r=a to r=0 of [ -((V/a) * r) * (8 * pi * epsilon_0 * (V/a) * r * dr) ]")
    print("    H = integral from r=a to r=0 of [ -8 * pi * epsilon_0 * (V**2/a**2) * r**2 * dr ]\n")

    # Step 4: Solve the integral and state the final answer
    print("--- Step 4: Final Result ---")
    print("Evaluating the integral:")
    print("    H = -8 * pi * epsilon_0 * (V**2/a**2) * [r**3 / 3]_from_a_to_0")
    print("    H = -8 * pi * epsilon_0 * (V**2/a**2) * (0**3 / 3 - a**3 / 3)")
    print("    H = -8 * pi * epsilon_0 * (V**2/a**2) * (-a**3 / 3)")
    print("    H = (8/3) * pi * epsilon_0 * (V**2 * a**3) / a**2\n")
    print("Simplifying gives the final equation for the total Joule heat dissipated.")
    
    # Print the final equation with numbers
    numerator = 8
    denominator = 3
    print("The final equation is:")
    print(f"H = ({numerator}/{denominator}) * pi * epsilon_0 * a * V**2")


solve_joule_heat()
>>> H = (8/3) * pi * epsilon_0 * a * V**2