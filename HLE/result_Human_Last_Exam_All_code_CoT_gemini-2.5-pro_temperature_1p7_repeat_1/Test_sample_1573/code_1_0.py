import numpy as np

def check_concyclic():
    """
    Checks if five given points are concyclic.
    A circle is defined by three points. The equation of a circle can be
    written as x^2 + y^2 + ax + by + c = 0. We will use three points to find
    the coefficients a, b, and c, and then check if the other two points
    satisfy the equation.
    """
    p1 = (0, 0)
    p2 = (2, 0)
    p3 = (2, 2)
    p4 = (0, 2)
    p5 = (1, 4)
    
    points = [p1, p2, p3, p4, p5]
    
    print(f"The five base points of the chair legs are: {points}\n")
    
    # We will use three non-collinear points to define the circle.
    # Let's use P1, P2, and P5.
    # P1(0,0), P2(2,0), P5(1,4)
    
    # For a point (x, y) on the circle, we have:
    # a*x + b*y + c = -(x^2 + y^2)
    # This gives a system of linear equations for a, b, c.
    
    # Using P1=(0,0): a*0 + b*0 + c = -(0^2+0^2)  => c = 0
    # Using P2=(2,0): a*2 + b*0 + c = -(2^2+0^2)  => 2a + c = -4
    # Using P5=(1,4): a*1 + b*4 + c = -(1^2+4^2)  => a + 4b + c = -17
    
    # From the first equation, c = 0.
    # Substitute c=0 into the second: 2a = -4 => a = -2
    # Substitute a=-2 and c=0 into the third: -2 + 4b = -17 => 4b = -15 => b = -3.75
    
    c = 0.0
    a = -2.0
    b = -15.0 / 4.0
    
    print("Step 1: Determine the equation of the circle passing through the first, second, and fifth points.")
    print(f"The first point P1=(0,0) gives the equation: a*0 + b*0 + c = 0. This implies c = 0.")
    print(f"The second point P2=(2,0) gives the equation: a*2 + b*0 + 0 = -(2^2 + 0^2). This implies 2a = -4, so a = {a}.")
    print(f"The fifth point P5=(1,4) gives the equation: a*1 + b*4 + 0 = -(1^2 + 4^2). This implies a + 4b = -17.")
    print(f"Substituting a = {a}, we get {a} + 4b = -17, which means 4b = -15, so b = {b}.")
    
    print(f"\nThe equation of the circle is: x^2 + y^2 + ({a})x + ({b})y + ({c}) = 0")
    print(f"Which simplifies to: x^2 + y^2 - 2x - 3.75y = 0\n")

    print("Step 2: Check if the remaining points lie on this circle.")
    
    # Check P3 = (2,2)
    x, y = p3
    result_p3 = x**2 + y**2 + a*x + b*y + c
    print(f"Testing point P3={p3}:")
    print(f" plugging into the equation: ({x})^2 + ({y})^2 + ({a})*({x}) + ({b})*({y}) + ({c}) = {result_p3}")
    if np.isclose(result_p3, 0):
        print("Point P3 is on the circle.")
    else:
        print(f"Point P3 is NOT on the circle, since the result is {result_p3}, not 0.")

    # Check P4 = (0,2)
    x, y = p4
    result_p4 = x**2 + y**2 + a*x + b*y + c
    print(f"\nTesting point P4={p4}:")
    print(f" plugging into the equation: ({x})^2 + ({y})^2 + ({a})*({x}) + ({b})*({y}) + ({c}) = {result_p4}")
    if np.isclose(result_p4, 0):
        print("Point P4 is on the circle.")
    else:
        print(f"Point P4 is NOT on the circle, since the result is {result_p4}, not 0.")

    print("\nConclusion: Since at least one point does not lie on the circle defined by the other three, the five points are not concyclic.")

check_concyclic()