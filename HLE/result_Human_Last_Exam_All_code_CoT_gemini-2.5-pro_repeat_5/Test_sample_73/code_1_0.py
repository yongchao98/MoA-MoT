import sympy as sp

def solve_genus():
    """
    Calculates the genus of the configuration space of a hinged pentagon
    with two adjacent vertices fixed.

    The method uses Morse theory on the parameter space of the linkage.
    1. The configuration is parameterized by two angles, x and y, on a torus.
    2. A function f(x, y) is defined for the squared distance between the
       key vertices, which determines the allowed configuration region M.
    3. The genus is calculated as g = 1 - chi(M), where chi(M) is the
       Euler characteristic of M.
    4. chi(M) is found by finding the critical points of f and summing their
       Morse indices for those points inside M.
    """
    print("Step 1: Define the distance-squared function f(x, y) on the parameter torus.")
    # x and y are the angles theta_3 and theta_5 from the derivation.
    # The function f is the squared distance between vertices V3 and V5, normalized by L^2.
    # f(x,y) = |V3-V5|^2 / L^2
    x, y = sp.symbols('x y')
    f = 3 + 2*sp.cos(x) - 2*sp.cos(y) - 2*sp.cos(x - y)
    print(f"f(x, y) = {f}\n")

    print("Step 2: Find the critical points of f(x, y) by solving grad(f) = 0.")
    # The critical points are found by solving the system:
    # sin(x) = sin(x-y)
    # sin(y) = sin(x-y)
    # This leads to sin(x) = sin(y), which means x=y or x=pi-y.
    # Solving these cases gives the following 6 unique points on the torus.
    crit_points = [
        {'name': 'P1', 'coords': (sp.pi * 0, sp.pi * 0)},
        {'name': 'P2', 'coords': (sp.pi, sp.pi)},
        {'name': 'P3', 'coords': (sp.pi, sp.pi * 0)},
        {'name': 'P4', 'coords': (2*sp.pi/3, sp.pi/3)},
        {'name': 'P5', 'coords': (4*sp.pi/3, 5*sp.pi/3)},
        {'name': 'P6', 'coords': (sp.pi * 0, sp.pi)},
    ]
    print("Found 6 critical points.\n")

    print("Step 3: Classify critical points using the Hessian matrix and check if they are in region M (f <= 4).")
    # The region M corresponds to valid configurations, where the distance between
    # V3 and V5 is at most 2L. This means f(x,y) <= 4.
    
    # Compute Hessian determinant
    f_x = sp.diff(f, x)
    f_y = sp.diff(f, y)
    f_xx = sp.diff(f_x, x)
    f_xy = sp.diff(f_x, y)
    f_yy = sp.diff(f_y, y)
    H_det = f_xx * f_yy - f_xy**2

    num_minima_inside = 0
    num_saddles_inside = 0

    print("-" * 65)
    print(f"{'Point':<5} | {'f value':<10} | {'Type':<10} | {'Index':<6}| {'In M (f<=4)?'}")
    print("-" * 65)
    
    for cp in crit_points:
        px, py = cp['coords']
        f_val = f.subs({x: px, y: py})
        det_val = H_det.subs({x: px, y: py})
        fxx_val = f_xx.subs({x: px, y: py})

        type_str = ""
        index = None
        if det_val > 0:
            if fxx_val > 0:
                type_str = "Minimum"
                index = 0
            else:
                type_str = "Maximum"
                index = 2
        elif det_val < 0:
            type_str = "Saddle"
            index = 1
        else:
            type_str = "Degenerate"
            index = 'N/A'

        is_inside = "Yes" if f_val <= 4 else "No"
        
        if is_inside == "Yes":
            if index == 0:
                num_minima_inside += 1
            elif index == 1:
                num_saddles_inside += 1
        
        print(f"{cp['name']:<5} | {f_val:<10.2f} | {type_str:<10} | {str(index):<6}| {is_inside}")
    print("-" * 65)

    print("\nStep 4: Calculate the Euler characteristic chi(M).")
    chi_M = num_minima_inside - num_saddles_inside
    print(f"The critical points inside M consist of {num_minima_inside} minima (index 0) and {num_saddles_inside} saddles (index 1).")
    print(f"chi(M) = (Number of minima) - (Number of saddles)")
    print(f"chi(M) = {num_minima_inside} - {num_saddles_inside} = {chi_M}\n")

    print("Step 5: Calculate the genus g of the configuration space.")
    genus = 1 - chi_M
    print("The genus g is related to chi(M) by the formula g = 1 - chi(M).")
    print(f"g = 1 - ({chi_M}) = {genus}")

    print("\nFinal Answer:")
    print(f"<<<{genus}>>>")

if __name__ == '__main__':
    solve_genus()