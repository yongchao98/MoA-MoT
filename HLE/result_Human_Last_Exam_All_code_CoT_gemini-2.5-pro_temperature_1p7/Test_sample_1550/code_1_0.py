import math
from fractions import Fraction

def solve_intersection():
    """
    Calculates the number of real intersection points between two conics.
    C1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    C2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    """
    # Coefficients of the conics C1 and C2 as Fractions
    C1_coeffs = {'A': Fraction(164), 'B': Fraction(-216), 'C': Fraction(72),
                 'D': Fraction(-16), 'E': Fraction(0), 'F': Fraction(31)}
    C2_coeffs = {'A': Fraction(864), 'B': Fraction(-1056), 'C': Fraction(324),
                 'D': Fraction(-560), 'E': Fraction(324), 'F': Fraction(149)}

    print("Step 1: Find a k that makes C1 + k*C2 a parabola.")
    print("A conic Ax^2 + Bxy + Cy^2 + ... = 0 is a parabola if B^2 - 4AC = 0.")
    print("For the pencil C1 + k*C2, the coefficients are:")
    print("A(k) = 164 + 864k")
    print("B(k) = -216 - 1056k")
    print("C(k) = 72 + 324k")
    print("\nThe condition B(k)^2 - 4A(k)C(k) = 0 gives a quadratic equation in k.")
    
    # Symbolic expansion of B(k)^2 - 4A(k)C(k)
    # (-216 - 1056k)^2 - 4(164 + 864k)(72 + 324k)
    # (46656 + 456192k + 1115136k^2) - 4(11808 + 115344k + 279936k^2)
    # 46656 - 47232 + (456192 - 461376)k + (1115136 - 1119744)k^2
    # -576 - 5184k - 4608k^2 = 0
    # Dividing by -576 gives: 1 + 9k + 8k^2 = 0
    
    k_eq_a = Fraction(8)
    k_eq_b = Fraction(9)
    k_eq_c = Fraction(1)
    
    print(f"The equation for k is: {int(k_eq_a)}*k^2 + {int(k_eq_b)}*k + {int(k_eq_c)} = 0")
    print("Factoring gives: (8k + 1)(k + 1) = 0")
    
    k1 = Fraction(-1, 8)
    k2 = Fraction(-1)
    print(f"The solutions for k are k = {k1} and k = {k2}.\n")
    
    print(f"Step 2: Choose k = {k1} and form the conic 8*C1 - C2 = 0.")
    k = k1
    
    # 8*C1 - C2
    coeffs = {key: 8 * C1_coeffs[key] - C2_coeffs[key] for key in C1_coeffs}
    
    print("The equation of the resulting parabola P is:")
    print(f"{int(coeffs['A'])}*x^2 + {int(coeffs['B'])}*x*y + {int(coeffs['C'])}*y^2 + {int(coeffs['D'])}*x + {int(coeffs['E'])}*y + {int(coeffs['F'])} = 0")
    
    print("\nThe quadratic part is 448x^2 - 672xy + 252y^2 = 28 * (16x^2 - 24xy + 9y^2) = 28 * (4x - 3y)^2.")
    print("The linear part is 432x - 324y = 108 * (4x - 3y).")
    print("Let u = 4x - 3y. The equation for P becomes a quadratic in u:")
    
    # 28u^2 + 108u + 99 = 0
    u_eq_a = Fraction(28)
    u_eq_b = Fraction(108)
    u_eq_c = Fraction(99)
    
    print(f"{int(u_eq_a)}*u^2 + {int(u_eq_b)}*u + {int(u_eq_c)} = 0\n")
    
    print("Step 3: Solve for u.")
    u_discriminant = u_eq_b**2 - 4 * u_eq_a * u_eq_c
    u_sqrt_d = Fraction(math.isqrt(u_discriminant.numerator), math.isqrt(u_discriminant.denominator))
    
    u1 = (-u_eq_b + u_sqrt_d) / (2 * u_eq_a)
    u2 = (-u_eq_b - u_sqrt_d) / (2 * u_eq_a)

    print(f"The solutions are u1 = {u1} and u2 = {u2}.")
    print("This means the intersection points must lie on two parallel lines:")
    print(f"L1: 4x - 3y = {u1}  (or 8x - 6y + 3 = 0)")
    print(f"L2: 4x - 3y = {u2}  (or 56x - 42y + 33 = 0)\n")
    
    total_intersections = 0
    
    # We substitute y = (4x - u)/3 into C1
    # 164x^2 - 216x(4x-u)/3 + 72((4x-u)/3)^2 - 16x + 31 = 0
    # ... simplifies to ...
    # 4x^2 - 8(u+2)x + (8u^2 + 31) = 0

    lines = [
        {'name': 'L1', 'u': u1},
        {'name': 'L2', 'u': u2}
    ]
    
    for line in lines:
        u = line['u']
        name = line['name']
        print(f"Step 4: Intersect C1 with line {name}.")
        
        # Coefficients of the resulting quadratic equation in x: ax^2 + bx + c = 0
        a = Fraction(4)
        b = -8 * (u + 2)
        c = 8 * u**2 + 31
        
        print(f"Substituting y from {name} into C1 gives the quadratic equation:")
        print(f"({a})*x^2 + ({b})*x + ({c}) = 0")
        
        discriminant = b**2 - 4 * a * c
        
        print(f"The discriminant is D = {discriminant}.")
        
        if discriminant > 0:
            intersections = 2
            print(f"D > 0, so there are 2 real intersection points on this line.\n")
        elif discriminant == 0:
            intersections = 1
            print(f"D = 0, so the line is tangent. There is 1 real intersection point on this line.\n")
        else:
            intersections = 0
            print(f"D < 0, so there are 0 real intersection points on this line.\n")
            
        total_intersections += intersections
        
    print("--------------------------------------------------")
    print(f"Total number of real intersection points is {total_intersections}.")

solve_intersection()