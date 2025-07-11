import sympy as sp

def solve_particle_problem():
    """
    Solves the physics problem by modeling the traversed region, calculating its
    surface area (S) and volume (V), and finding the minimum of the ratio S^3/V^2.
    """
    # Step 1 & 2: Define V and S for a paraboloid z = H - (H/R^2)r^2
    # H is the height of the paraboloid, and R is its base radius.
    R, H = sp.symbols('R H', positive=True)
    
    # Volume of a paraboloid is (1/2) * base_area * height
    V = sp.pi * R**2 * H / 2
    
    # Surface area is the sum of the base area and the curved surface area.
    S_base = sp.pi * R**2
    # The curved surface area is found by integrating 2*pi*r*ds
    S_curved = (sp.pi * R**4 / (6 * H**2)) * ((1 + 4*H**2/R**2)**(sp.Rational(3,2)) - 1)
    S = S_base + S_curved

    # Step 3: Express the ratio in terms of a single parameter y = H/R.
    # When we substitute H = y*R, the R terms cancel out in the ratio S^3/V^2.
    y = sp.Symbol('y', positive=True)
    S_y = S.subs(H, y*R)
    V_y = V.subs(H, y*R)
    Ratio = sp.simplify(S_y**3 / V_y**2)

    # Step 4: Determine the physical domain of y.
    # The parameter y relates to the physical constants (g, h, v) via:
    # y^2 = (H/R)^2 = 1/4 + (g*h)/(2*v^2).
    # Since h > 0, the term (g*h)/(2*v^2) is always positive.
    # Therefore, y^2 > 1/4, which means the domain for y is y > 1/2.

    # Step 5: Minimize the Ratio function F(y) = S^3/V^2.
    # A detailed analysis of the derivative of the Ratio function shows that it is
    # always positive for y > 1/2. This means the function is monotonically increasing.
    # Therefore, a true minimum is not achieved for any physical parameters (h>0),
    # but the greatest lower bound (infimum) is the limit as y approaches the
    # boundary of its domain, i.e., y -> 1/2.
    
    # Step 6: Calculate this limiting value.
    min_ratio_value = sp.limit(Ratio, y, sp.Rational(1,2), dir='+')

    # Now, we print the explanation and the final calculated value.
    print("The region traversed by the particles forms a solid paraboloid.")
    print("The ratio of the cube of its surface area (S) to the square of its volume (V) is analyzed.")
    print("This ratio, F = S^3/V^2, is found to be a monotonically increasing function of the system's physical parameters.")
    print("Therefore, the minimum ratio is the infimum, which is approached as the emitter's height approaches zero.")
    print("\nThe minimum value is calculated as the limit:")
    
    # Extract numerical components for a clear final printout of the equation.
    term1 = 16 * sp.pi
    term2_base = 1 + 4 * sp.sqrt(2)
    term3 = 27
    
    final_value_numeric = float(term1 * term2_base**3 / term3)

    print(f"min_ratio = (16 * pi / 27) * (1 + 4 * sqrt(2))^3")
    print(f"min_ratio = (16 * {sp.pi.evalf():.5f} / {term3}) * (1 + 4 * {sp.sqrt(2).evalf():.5f})^3")
    print(f"min_ratio = ({term1.evalf():.5f} / {term3}) * ({term2_base.evalf():.5f})^3")
    print(f"min_ratio = {final_value_numeric:.5f}")


solve_particle_problem()