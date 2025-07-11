import math

def solve_overhang():
    """
    Calculates the integers a, b, c for the maximal overhang problem.
    """
    # Based on the plan, we solve for the CM positions.
    # Let d = 1/sqrt(2) be the maximal horizontal displacement.
    # x1 = x2 + d
    # (x1 + x2)/2 = x3 + d  => x2 + d/2 = x3 + d => x2 = x3 + d/2
    # (x1 + x2 + x3)/3 = 0 => ((x3 + 3d/2) + (x3 + d/2) + x3)/3 = 0
    # => (3*x3 + 2d)/3 = 0 => x3 + 2d/3 = 0 => x3 = -2d/3
    
    d = 1/math.sqrt(2)
    
    x3 = -2 * d / 3
    x2 = x3 + d / 2
    x1 = x2 + d
    
    # The maximal overhang is achieved by also rotating the top cube C1.
    # The rightmost point of a rotated cube i is at xi + d.
    overhang1 = x1 + d
    overhang2 = x2 + d
    overhang3 = x3 + d
    
    max_overhang = max(overhang1, overhang2, overhang3)
    
    # The maximal overhang is 11*d/6 = 11/(6*sqrt(2)) = 11*sqrt(2)/12.
    # We need to write this as (a + sqrt(b))/(1+c).
    # overhang = (11 * sqrt(2)) / 12 = sqrt(121 * 2) / 12 = sqrt(242) / 12
    # This is (0 + sqrt(242)) / 12.
    
    a = 0
    b = 242  # b = 11^2 * 2
    # Denominator is 12, so 1 + c = 12
    c = 11
    
    print("The maximal overhang is achieved by rotating all three cubes by 45 degrees.")
    print(f"The calculated maximal overhang is (11 * sqrt(2)) / 12.")
    print(f"This can be written in the specified format (a + sqrt(b))/(1 + c) as:")
    print(f"({a} + sqrt({b}))/(1 + {c})")
    print("\nThe values for a, b, and c are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")

solve_overhang()