import math

def solve_triangle_median_range():
    """
    Calculates the range of values for the median m given the height AD,
    angle bisector AE, and the condition that angle A is acute.
    """

    # Given values from the problem statement
    # AD is the height on side BC
    h_AD = 12
    # AE is the angle bisector of angle A
    l_AE = 13

    # We can model this problem in a 2D coordinate system. Let D be the origin (0,0).
    # A is at (0, 12). E and F are on the x-axis.
    # Triangle ADE is a right-angled triangle with the right angle at D.
    # We can find the length of the side DE using the Pythagorean theorem.
    # DE^2 = AE^2 - AD^2
    DE_sq = l_AE**2 - h_AD**2
    
    # DF is the distance from the foot of the altitude D to the foot of the median F.
    # AF is the median, with length m. Triangle ADF is also a right-angled triangle.
    # m^2 = AD^2 + DF^2
    # DF^2 = m^2 - AD^2

    # The condition that angle A is acute, combined with the other geometric constraints,
    # implies that the foot of the median (F) must lie between the foot of the
    # altitude (D) and the foot of the angle bisector (E).
    # This means the distance DF must be greater than 0 and less than the distance DE.
    # 0 < DF^2 < DE^2

    # Substitute the expressions for DF^2 and DE^2 into the inequality:
    # 0 < m^2 - AD^2 < AE^2 - AD^2
    
    # Let's plug in the numbers:
    # 0 < m^2 - 12^2 < 13^2 - 12^2
    # 0 < m^2 - 144 < 169 - 144
    # 0 < m^2 - 144 < 25

    # To find the range for m, we first solve for m^2 by adding 144 to all parts:
    # 144 < m^2 < 25 + 144
    # 144 < m^2 < 169

    # Now, we take the square root of all parts. Since m is a length, it must be positive.
    lower_bound_m_sq = 144
    upper_bound_m_sq = 169
    
    lower_bound_m = math.sqrt(lower_bound_m_sq)
    upper_bound_m = math.sqrt(upper_bound_m_sq)

    # The final equation is the inequality describing the range of m.
    # We print each number in this final equation.
    print(f"The condition that angle A is acute requires the median length m to be in the range:")
    print(f"{int(lower_bound_m)} < m < {int(upper_bound_m)}")

solve_triangle_median_range()