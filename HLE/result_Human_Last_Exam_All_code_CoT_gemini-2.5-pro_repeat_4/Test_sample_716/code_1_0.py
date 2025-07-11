import math

def solve_hausdorff_dimension():
    """
    This script explains the analytical derivation of the Hausdorff dimension
    for the given curve.
    """
    print("Problem: Find the Hausdorff dimension of the curve parametrized by:")
    print("x(t) = sin(pi * t)")
    print("y(t) = sin(t)")
    print("z(t) = cos(2t)")
    print("-" * 30)

    print("Step 1: Analyze the smoothness of the curve.")
    print("A curve is smooth if its derivative vector r'(t) is continuous and never the zero vector.")
    print("\nThe vector parametrization is r(t) = (sin(pi*t), sin(t), cos(2t)).")
    print("-" * 30)

    print("Step 2: Calculate the derivative vector r'(t).")
    print("The derivative of each component with respect to t is:")
    print("x'(t) = d/dt(sin(pi*t)) = pi * cos(pi*t)")
    print("y'(t) = d/dt(sin(t))     = cos(t)")
    print("z'(t) = d/dt(cos(2t))    = -2 * sin(2t)")
    print("\nSo, r'(t) = (pi * cos(pi*t), cos(t), -2 * sin(2t)).")
    print("The components are all continuous functions, so r'(t) is continuous.")
    print("-" * 30)

    print("Step 3: Check if r'(t) can be the zero vector (0, 0, 0).")
    print("This requires all three components to be zero for the same value of t.")
    print("1) pi * cos(pi*t) = 0  =>  pi*t = pi/2 + n*pi  =>  t = 1/2 + n (for integer n)")
    print("2) cos(t) = 0        =>  t = pi/2 + k*pi      (for integer k)")
    print("3) -2 * sin(2t) = 0    =>  2t = m*pi          =>  t = m*pi/2       (for integer m)")
    print("\nLet's see if condition 1 and 2 can be met simultaneously:")
    print("1/2 + n = pi/2 + k*pi")
    print("1/2 + n = (1/2 + k) * pi")
    print("\nThe left side of this equation is a rational number. The right side is an irrational number (a rational number times pi).")
    print("A rational number cannot equal an irrational number. Therefore, there is no value of t for which both x'(t) and y'(t) are zero.")
    print("\nSince it's impossible for even two components to be zero simultaneously, the derivative vector r'(t) is never the zero vector.")
    print("-" * 30)

    print("Step 4: Conclusion.")
    print("The curve is smooth because its derivative is continuous and never zero.")
    print("A fundamental theorem of fractal geometry states that any smooth, rectifiable curve has a Hausdorff dimension of 1.")
    print("\nTherefore, the Hausdorff dimension of the given curve is 1.")
    print("-" * 30)
    
    dimension = 1
    print(f"The final equation is: Hausdorff Dimension = {dimension}")

if __name__ == "__main__":
    solve_hausdorff_dimension()
