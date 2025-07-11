import math

def solve_toilet_paper_problem():
    """
    Calculates the time it takes for a falling toilet paper roll to completely unroll.
    """
    # --- Given Parameters ---
    d_cylinder_cm = 4.0   # cm
    thickness_mm = 0.5    # mm
    N_wraps = 100         # number of wraps
    M_paper_g = 200.0     # grams
    M_cylinder_g = 20.0   # grams
    g = 9.8               # m/s^2

    # --- Convert to SI Units ---
    d_cylinder = d_cylinder_cm / 100.0  # m
    thickness = thickness_mm / 1000.0   # m
    M_paper = M_paper_g / 1000.0        # kg
    M_cylinder = M_cylinder_g / 1000.0  # kg

    # --- Step 1: Calculate derived parameters ---
    r_cylinder = d_cylinder / 2
    M_total = M_paper + M_cylinder

    # --- Step 2: Calculate the total length of the paper (L) ---
    # The length is the sum of the circumferences of N concentric circles.
    # The sum of radii is an arithmetic series: N*r_cylinder + thickness * N*(N+1)/2
    sum_of_radii = N_wraps * r_cylinder + thickness * (N_wraps * (N_wraps + 1) / 2)
    L = 2 * math.pi * sum_of_radii

    # --- Step 3: Calculate the time using the analytical solution ---
    # The analytical solution for the total time 't' is:
    # t = sqrt( (M_total * L) / (2 * g * M_paper) ) * ( asin((M_paper - M_cylinder) / M_total) + pi/2 )

    # Term 1: The coefficient in front of the angle term
    coeff = math.sqrt((M_total * L) / (2 * g * M_paper))

    # Term 2: The argument of the arcsin function
    asin_arg = (M_paper - M_cylinder) / M_total
    
    # Term 3: The angle part of the formula
    angle_term = math.asin(asin_arg) + math.pi / 2

    # Final time calculation
    time = coeff * angle_term

    # --- Step 4: Print the detailed calculation and final answer ---
    print("--- Toilet Paper Unrolling Calculation ---")
    print(f"\n1. Calculated Parameters:")
    print(f"   - Total Paper Length (L): {L:.4f} m")
    print(f"   - Total Initial Mass (M_total): {M_total:.3f} kg")

    print("\n2. Final Time Equation:")
    print("   t = sqrt( (M_total * L) / (2 * g * M_paper) ) * [arcsin((M_paper - M_cylinder) / M_total) + π/2]")
    
    print("\n3. Plugging in the numbers:")
    print(f"   t = sqrt( ({M_total:.3f} kg * {L:.4f} m) / (2 * {g} m/s^2 * {M_paper:.3f} kg) ) * [arcsin(({M_paper:.3f} - {M_cylinder:.3f}) / {M_total:.3f}) + π/2]")
    print(f"   t = sqrt({coeff**2:.4f}) * [arcsin({asin_arg:.4f}) + {math.pi/2:.4f}]")
    print(f"   t = {coeff:.4f} * [{math.asin(asin_arg):.4f} + {math.pi/2:.4f}]")
    print(f"   t = {coeff:.4f} * {angle_term:.4f}")
    
    print("\n--- Final Answer ---")
    print(f"The time it takes for the toilet paper to fully unroll is {time:.2f} seconds.")

# Run the solver
solve_toilet_paper_problem()
print("\n<<<3.19>>>")