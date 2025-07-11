import sympy

def check_smoothness(func, var, point, max_order=5):
    """
    Checks the differentiability of a function at a point up to a certain order.
    """
    print(f"\n--- Analyzing smoothness of '{func}' at {var}={point} ---")
    is_smooth = True
    for k in range(1, max_order + 1):
        try:
            # Calculate the k-th derivative
            deriv_expr = sympy.diff(func, var, k)
            
            # Evaluate the derivative at the point. This might fail.
            # A more robust way is to check the limits from both sides.
            lim_pos = sympy.limit(deriv_expr, var, point, dir='+')
            lim_neg = sympy.limit(deriv_expr, var, point, dir='-')

            print(f"Order {k} derivative expression: {deriv_expr}")
            
            if lim_pos.is_real and lim_neg.is_real and lim_pos == lim_neg:
                print(f"Order {k} derivative at {point} is: {lim_pos}")
            else:
                print(f"Order {k} derivative does not exist at {point}.")
                print(f"  Limit from above: {lim_pos}")
                print(f"  Limit from below: {lim_neg}")
                is_smooth = False
                break
        except Exception as e:
            print(f"Could not compute order {k} derivative at {point}. Error: {e}")
            is_smooth = False
            break
    if is_smooth:
        print(f"The function appears to be smooth at {point} (checked up to order {max_order}).")

# --- Main execution ---
x = sympy.Symbol('x')

# The function |x| is not smooth because its first derivative is not defined at x=0.
f1 = sympy.Abs(x)
check_smoothness(f1, x, 0)

# The function |x^2| = x^2 is smooth.
f2 = sympy.Abs(x**2)
check_smoothness(f2, x, 0)

# The function |x^3| is C^2 but not C^3 at x=0. Its 3rd derivative does not exist.
# This illustrates why a parameterization like (t^3, |t^3|) would not be a smooth curve.
f3 = sympy.Abs(x**3)
check_smoothness(f3, x, 0)