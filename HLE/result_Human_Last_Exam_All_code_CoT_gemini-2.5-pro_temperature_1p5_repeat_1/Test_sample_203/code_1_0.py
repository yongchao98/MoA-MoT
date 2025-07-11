import math

def solve_triangle_median_range():
    """
    Calculates the range for the median m based on the derived formulas.
    """
    # Given values
    h_a = 12 # Altitude AD
    t_a = 13 # Angle bisector AE
    
    # Derived values
    # DE = sqrt(AE^2 - AD^2)
    k = math.sqrt(t_a**2 - h_a**2) # k = 5
    
    # Derived relations for d = x2 - x1
    # 119*d = 1440 + 10*P, where P = x1*x2
    # For angle A to be acute, 0 < P < 144
    
    # Lower bound for d (as P -> 0)
    d_lower_num = 1440
    d_den = 119
    
    # Upper bound for d (as P -> 144)
    d_upper_num = 2880
    
    print("Step 1: The condition for angle A to be acute translates to a range for a quantity P, where 0 < P < 144.")
    print("Step 2: The geometric setup gives a relation between another quantity d and P: 119*d = 1440 + 10*P.")
    print(f"Step 3: This implies a range for d: ({d_lower_num}/{d_den}) < d < ({d_upper_num}/{d_den}).")
    
    # Range for the median m, where m^2 = (d/2)^2 + h_a^2
    
    # Calculate Lower bound for m
    # m_lower^2 = (d_lower/2)^2 + h_a^2 = (720/119)^2 + 144
    m_lower_sq_num = 720**2 + h_a**2 * d_den**2
    # m_lower_sq_num = 144 * 17761
    m_lower_sqrt_arg = 17761
    m_lower_num = h_a * math.sqrt(m_lower_sqrt_arg)
    m_lower_den = d_den
    m_lower_val = m_lower_num / m_lower_den
    
    # Calculate Upper bound for m
    # m_upper^2 = (d_upper/2)^2 + h_a^2 = (1440/119)^2 + 144
    m_upper_num = 12 * 169
    m_upper_den = 119
    m_upper_val = m_upper_num / m_upper_den
    
    print(f"Step 4: The median m is related to d by m^2 = (d/2)^2 + {h_a**2}.")
    print("This gives the final range for m.")
    
    print("\nThe range for m is:")
    print(f"(12 * sqrt({m_lower_sqrt_arg})) / {m_lower_den} < m < {m_upper_num} / {m_upper_den}")
    print("\nNumerically, this is:")
    print(f"{m_lower_val:.4f} < m < {m_upper_val:.4f}")

solve_triangle_median_range()