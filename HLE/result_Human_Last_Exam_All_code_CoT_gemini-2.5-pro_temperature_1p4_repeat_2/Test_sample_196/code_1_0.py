import sympy as sp

def solve_volume_problem():
    """
    This function calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.
    
    The steps are:
    1. Define the parameters of the cone and ellipsoid.
    2. Determine the plane of tangency and the radius of the circular intersection.
    3. Calculate the volume of the cone tip above this plane.
    4. Calculate the volume of the ellipsoid cap above this plane using integration.
    5. The final enclosed volume is the difference between the cone volume and the ellipsoid volume.
    """
    
    # --- Step 1 & 2: Define parameters and geometry ---
    # Cone vertex is at (0, 4, 0).
    # Ellipsoid is x^2/3 + y^2/4 + z^2/3 = 1.
    # The plane of tangency between the cone and ellipsoid is y = 1.
    # The intersection is a circle x^2 + z^2 = R^2 on the plane y = 1.
    # Substituting y=1 into the ellipsoid equation: x^2/3 + 1/4 + z^2/3 = 1 => x^2/3 + z^2/3 = 3/4 => x^2 + z^2 = 9/4.
    R_base_sq = sp.Rational(9, 4)
    R_base = sp.sqrt(R_base_sq)

    # --- Step 3: Calculate the volume of the cone tip ---
    # The cone tip has its vertex at y=4 and its base on the plane y=1.
    h_cone = 4 - 1
    # Volume of a cone is (1/3) * pi * R^2 * h
    V_cone_tip = sp.Rational(1, 3) * sp.pi * R_base_sq * h_cone

    # --- Step 4: Calculate the volume of the ellipsoid cap ---
    # The cap is the part of the ellipsoid from y=1 to its top at y=2.
    # We integrate the cross-sectional area A(y) from y=1 to y=2.
    # The cross-section at height y is a circle with squared radius r(y)^2 = 3 * (1 - y^2/4).
    y = sp.Symbol('y')
    Area_y = 3 * sp.pi * (1 - y**2 / 4)
    # Integrate A(y) from y=1 to y=2 to get the volume of the cap.
    V_ellipsoid_cap = sp.integrate(Area_y, (y, 1, 2))

    # --- Step 5: Calculate the final enclosed volume ---
    V_enclosed = V_cone_tip - V_ellipsoid_cap

    # --- Output the results step-by-step ---
    print("--- Calculation of the Enclosed Volume ---")
    print("\nThe enclosed volume is the volume of the cone tip minus the volume of the ellipsoid cap.")
    
    print("\nStep 1: Calculate the Volume of the Cone Tip (V_cone)")
    print(f"The cone has height h = 4 - 1 = {h_cone}")
    print(f"The base radius is R = {R_base}")
    print(f"V_cone = (1/3) * pi * ({R_base})^2 * ({h_cone}) = {V_cone_tip}")

    print("\nStep 2: Calculate the Volume of the Ellipsoid Cap (V_ellipsoid)")
    print(f"The volume is the integral of the cross-sectional area A(y) from y=1 to y=2.")
    print(f"V_ellipsoid = integral( 3*pi*(1 - y^2/4), (y, 1, 2) ) = {V_ellipsoid_cap}")

    print("\nStep 3: Final Volume Calculation")
    # Extract the numerical coefficients for the final equation display
    v_cone_num = V_cone_tip / sp.pi
    v_ell_num = V_ellipsoid_cap / sp.pi
    v_final_num = V_enclosed / sp.pi
    
    print("Final Volume = V_cone - V_ellipsoid")
    print(f"Final Volume = ({v_cone_num})*pi - ({v_ell_num})*pi = ({v_final_num})*pi")
    print(f"\nThe exact enclosed volume is: {V_enclosed}")
    
solve_volume_problem()