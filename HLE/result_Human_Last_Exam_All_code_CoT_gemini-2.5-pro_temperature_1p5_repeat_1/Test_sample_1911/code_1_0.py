import sympy

def check_boundedness():
    """
    This function demonstrates the topological difference between L and S^n.
    L is shown to be unbounded, while S^n is bounded. Since a diffeomorphism
    preserves boundedness (a property related to compactness), L and S^n
    cannot be diffeomorphic. This shows that statement C is false.
    """
    
    print("Analyzing Statement C: L can be given a smooth structure so that it is diffeomorphic to S^n for any n.")
    print("-" * 80)
    
    # --- Analysis of L = {(x, y) | y = |x|} ---
    print("1. Analyzing the set L.")
    t = sympy.Symbol('t', real=True)
    
    # A parameterization of L is (t, |t|)
    x_L = t
    y_L = sympy.Abs(t)
    
    # Calculate the squared norm (distance from origin) for a point on L
    norm_sq_L = x_L**2 + y_L**2
    norm_sq_L_simplified = sympy.simplify(norm_sq_L)
    
    # The norm is the square root of the squared norm
    norm_L = sympy.sqrt(norm_sq_L_simplified)
    
    print(f"A point on L can be parameterized by t as: ({x_L}, {y_L})")
    print(f"The squared norm of a point on L is ||(t, |t|)||² = {x_L}² + ({y_L})² = {norm_sq_L_simplified}")
    # We use a final equation with numbers to comply with the prompt's request.
    # The number is sqrt(2).
    print(f"The norm is ||(t, |t|)|| = {norm_L} = {sympy.sqrt(2)}*|t|.")
    print("As t approaches infinity, the norm also approaches infinity.")
    print("This shows that the set L is unbounded.\n")
    
    # --- Analysis of S^n (using S^1 as an example) ---
    print("2. Analyzing the set S^n (using n=1 as an example).")
    theta = sympy.Symbol('theta', real=True)
    
    # A parameterization of the unit circle S^1 is (cos(theta), sin(theta))
    x_S1 = sympy.cos(theta)
    y_S1 = sympy.sin(theta)
    
    # Calculate the squared norm for a point on S^1
    norm_sq_S1 = x_S1**2 + y_S1**2
    norm_sq_S1_simplified = sympy.simplify(norm_sq_S1)
    
    # The norm is the square root
    norm_S1 = sympy.sqrt(norm_sq_S1_simplified)
    
    print(f"A point on S^1 can be parameterized by theta as: ({x_S1}, {y_S1})")
    print(f"The squared norm of a point on S^1 is ||(cos(t), sin(t))||² = ({x_S1})² + ({y_S1})² = {norm_sq_S1_simplified}")
    print(f"The norm is ||(cos(t), sin(t))|| = {norm_S1}.")
    print("The norm is always 1, regardless of the value of theta.")
    print("This shows that the set S^1 is bounded.\n")
    
    # --- Conclusion ---
    print("3. Conclusion.")
    print("A diffeomorphism is a homeomorphism, which preserves topological properties like compactness.")
    print("Since L is unbounded (non-compact) and S^n is bounded (compact), they cannot be homeomorphic,")
    print("and therefore cannot be diffeomorphic.")
    print("Thus, statement C is false.")

if __name__ == '__main__':
    check_boundedness()
