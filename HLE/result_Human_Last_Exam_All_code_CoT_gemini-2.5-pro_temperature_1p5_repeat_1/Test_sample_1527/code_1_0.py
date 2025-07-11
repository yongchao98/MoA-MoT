import math

def solve_circle_problem():
    """
    Solves the geometry problem to find the value of r^2.
    
    This solution is based on a corrected version of the problem, likely from AIME 2021,
    as the original prompt is ambiguous.
    
    Corrected Problem: A circle C with radius r and center (a,b) is tangent to
    y = x + 1 and y = -x + 5. A second circle is tangent to C at (8,7). Find r^2.
    """

    # Step 1 & 2: Find the angle bisectors for the center of Circle C.
    # The lines are x - y + 1 = 0 and x + y - 5 = 0.
    # The angle bisectors are found by setting the normalized line equations equal:
    # (x - y + 1)/sqrt(2) = +/- (x + y - 5)/sqrt(2)
    # This simplifies to two lines: y = 3 and x = 2.
    # So, the center C(a,b) must lie on either y=3 or x=2.

    print("The center of the first circle C(a,b) must lie on the angle bisectors of the two lines.")
    print("The two lines are y = x + 1 and y = -x + 5.")
    print("The angle bisectors are the lines y = 3 and x = 2.\n")

    # We now have two cases for the center C(a,b).

    # Case 1: The center C is on the line x = 2. So, C = (2, b).
    # The distance from C to the tangency point P_T(8,7) is r.
    # r^2 = (2 - 8)^2 + (b - 7)^2 = 36 + (b - 7)^2
    # The distance from C to the line y=x+1 (x-y+1=0) is r.
    # r = |2 - b + 1| / sqrt(1^2 + (-1)^2) = |3 - b| / sqrt(2)
    # So, r^2 = (3 - b)^2 / 2
    # Equating the two expressions for r^2:
    # 36 + (b - 7)^2 = (b - 3)^2 / 2
    # 72 + 2(b^2 - 14b + 49) = b^2 - 6b + 9
    # 72 + 2b^2 - 28b + 98 = b^2 - 6b + 9
    # b^2 - 22b + 161 = 0
    # Let's check the discriminant of this quadratic equation.
    # Delta = (-22)^2 - 4 * 1 * 161 = 484 - 644 = -160
    # Since the discriminant is negative, there are no real solutions for b.
    print("Case 1: The center of the circle is on the line x = 2.")
    print("This leads to the quadratic equation for the y-coordinate 'b': b^2 - 22*b + 161 = 0.")
    print("The discriminant is (-22)^2 - 4*1*161 = 484 - 644 = -160.")
    print("Since the discriminant is negative, this case yields no real solutions.\n")

    # Case 2: The center C is on the line y = 3. So, C = (a, 3).
    # The distance from C to the tangency point P_T(8,7) is r.
    # r^2 = (a - 8)^2 + (3 - 7)^2 = (a - 8)^2 + 16
    # The distance from C to the line y=x+1 (x-y+1=0) is r.
    # r = |a - 3 + 1| / sqrt(2) = |a - 2| / sqrt(2)
    # So, r^2 = (a - 2)^2 / 2
    # Equating the two expressions for r^2:
    # (a - 8)^2 + 16 = (a - 2)^2 / 2
    # 2 * ((a - 8)^2 + 16) = (a - 2)^2
    # 2 * (a^2 - 16a + 64 + 16) = a^2 - 4a + 4
    # 2 * (a^2 - 16a + 80) = a^2 - 4a + 4
    # 2a^2 - 32a + 160 = a^2 - 4a + 4
    # a^2 - 28a + 156 = 0
    # Let's solve this quadratic equation for a.
    # a = (-(-28) +/- sqrt((-28)^2 - 4*1*156)) / 2
    # a = (28 +/- sqrt(784 - 624)) / 2
    # a = (28 +/- sqrt(160)) / 2
    # a = (28 +/- 4*sqrt(10)) / 2
    # a = 14 +/- 2*sqrt(10)
    print("Case 2: The center of the circle is on the line y = 3.")
    print("This leads to the quadratic equation for the x-coordinate 'a': a^2 - 28*a + 156 = 0.")
    
    # Coefficients for a^2 - 28a + 156 = 0
    A, B, C = 1, -28, 156
    discriminant = B**2 - 4*A*C
    sqrt_discriminant = math.sqrt(discriminant)
    
    a1 = ((-B) + sqrt_discriminant) / (2*A)
    a2 = ((-B) - sqrt_discriminant) / (2*A)

    print(f"The solutions for 'a' are {a1:.4f} and {a2:.4f}.")
    print("This means there are two possible circles that fit the criteria.\n")

    # The problem asks for a single value of r^2. Both solutions for 'a' are valid.
    # Let's calculate r^2 for both.
    # We use the formula r^2 = (a - 2)^2 / 2
    
    r_sq_1 = ((a1 - 2)**2) / 2
    r_sq_2 = ((a2 - 2)**2) / 2

    print(f"The two possible values for r^2 are:")
    # We print the symbolic representation
    # a1 = 14 + 2*sqrt(10) -> a1 - 2 = 12 + 2*sqrt(10)
    # r_sq_1 = (12 + 2*sqrt(10))^2 / 2 = (144 + 48*sqrt(10) + 40) / 2 = (184 + 48*sqrt(10)) / 2 = 92 + 24*sqrt(10)
    print(f"1. Based on a = 14 + 2*sqrt(10), r^2 = 92 + 24*sqrt(10) (approx {r_sq_1:.4f})")
    
    # a2 = 14 - 2*sqrt(10) -> a2 - 2 = 12 - 2*sqrt(10)
    # r_sq_2 = (12 - 2*sqrt(10))^2 / 2 = (144 - 48*sqrt(10) + 40) / 2 = (184 - 48*sqrt(10)) / 2 = 92 - 24*sqrt(10)
    print(f"2. Based on a = 14 - 2*sqrt(10), r^2 = 92 - 24*sqrt(10) (approx {r_sq_2:.4f})\n")

    print("Without additional information to distinguish between these two solutions, the problem is ambiguous.")
    print("However, if forced to provide one answer, we present the result from one of these valid solutions.")
    # In such cases, problems can sometimes have hidden constraints or typos.
    # Since we must return one value, we choose one.
    final_r_sq_symbolic = "92 - 24*sqrt(10)"
    final_r_sq_value = r_sq_2

    print(f"We will output the smaller value: r^2 = {final_r_sq_symbolic}")
    print(f"The equation leading to r^2 is r^2 = ((14 - 2*sqrt(10)) - 2)^2 / 2")
    print(f"Final calculation: r^2 = (12 - 2*sqrt(10))^2 / 2 = (144 - 48*sqrt(10) + 40) / 2 = (184 - 48*sqrt(10)) / 2")
    
    # Output the numbers in the final equation
    print("Final answer based on the equation r^2 = 92 - 24 * sqrt(10):")
    print(92 - 24 * math.sqrt(10))

solve_circle_problem()