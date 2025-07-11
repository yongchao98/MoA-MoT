import sympy

def find_shape_equation():
    """
    Analyzes the shape of the set S for a specific 2D example.
    We choose y1=[1,0] and y2=[1,1].
    The components of a vector in S are x1 and x2.
    We express them in terms of a parameter t (theta) and then eliminate t
    to find the implicit equation relating x1 and x2.
    """
    # Define symbols
    t = sympy.Symbol('t') # Represents theta
    x1 = sympy.Symbol('x1')
    x2 = sympy.Symbol('x2')

    # Expressions for x1 and x2 based on our example
    # y1 = [1, 0], y2 = [1, 1], s = [cos(t), sin(t)]
    # x1 = |<y1, s>|^2 = cos(t)^2
    # x2 = |<y2, s>|^2 = (cos(t) + sin(t))^2
    
    # We will derive the equation relating x1 and x2
    # From x1 = cos(t)**2, we get cos(t)**2 = x1 and sin(t)**2 = 1 - x1.
    # From x2 = (cos(t) + sin(t))**2 = 1 + 2*sin(t)*cos(t), we get 2*sin(t)*cos(t) = x2 - 1.
    # Squaring this: 4*sin(t)**2*cos(t)**2 = (x2 - 1)**2.
    # Substitute the expressions in terms of x1: 4*(1-x1)*x1 = (x2-1)**2
    
    # Let's form the polynomial equation from this derivation.
    # 4*x1 - 4*x1**2 = x2**2 - 2*x2 + 1
    # 4*x1**2 - 4*x1 + x2**2 - 2*x2 + 1 = 0
    
    equation = 4*x1**2 - 4*x1 + x2**2 - 2*x2 + 1
    
    # Extract coefficients of the quadratic equation A*x1^2 + B*x1*x2 + C*x2^2 + D*x1 + E*x2 + F = 0
    poly = sympy.poly(equation, x1, x2)
    
    A = poly.coeff_monomial(x1**2)
    B = poly.coeff_monomial(x1*x2)
    C = poly.coeff_monomial(x2**2)
    D = poly.coeff_monomial(x1)
    E = poly.coeff_monomial(x2)
    F = poly.coeff_monomial(1)
    
    print("For the example y1=[1,0], y2=[1,1], the shape is described by the equation:")
    print(f"{equation} = 0")
    print("\nThis is the equation of an ellipse, which is a 2D ellipsoid.")
    print("The general form is A*x1^2 + B*x1*x2 + C*x2^2 + D*x1 + E*x2 + F = 0.")
    print("The coefficients are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print(f"F = {F}")

if __name__ == '__main__':
    find_shape_equation()
