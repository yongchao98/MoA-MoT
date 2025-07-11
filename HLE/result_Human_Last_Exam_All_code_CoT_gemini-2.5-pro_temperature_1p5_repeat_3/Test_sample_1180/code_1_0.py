import sympy

def solve():
    """
    This problem asks for the thickness of a double point in the stable reduction of a hyperelliptic curve.
    
    The curve is C: z^2 = 2*x^5 + 2*x^3 + 1.
    The valuation v is normalized by v(2) = 1.
    
    1. Transform the curve to handle the branch point at infinity.
    Let x = 1/y. The new equation for the branch points becomes h(y) = y*(y^5 + 2*y^2 + 2) = 0.
    The branch points are now the roots of h(y).
    
    2. Partition the branch points.
    h(y) is factored over the integers (and thus over the 2-adic integers) into h1(y) = y and h2(y) = y^5 + 2*y^2 + 2.
    This factorization corresponds to a partition of the branch points into {infinity} and the set of 5 finite branch points.
    This partition gives rise to a stable model with two components connected by a node.
    
    3. Calculate the thickness.
    The thickness 't' of the node is the valuation of the resultant of the two polynomials in the factorization.
    t = v(Res(h1, h2)).
    
    Let's calculate the resultant of h1 and h2.
    h1(y) = y
    h2(y) = y^5 + 2*y^2 + 2
    
    The resultant of y and a polynomial P(y) is P(0).
    So, Res(y, y^5 + 2*y^2 + 2) = 0^5 + 2*0^2 + 2 = 2.
    """
    
    y = sympy.Symbol('y')
    h1 = y
    h2 = y**5 + 2*y**2 + 2
    
    # The resultant of a monomial y and a polynomial P(y) is P(0).
    resultant_value = h2.subs(y, 0)
    
    # The thickness is the valuation v(resultant), where v(2)=1.
    thickness = 1

    print(f"The equation for the branch points after transforming x=1/y is w^2 = y * (y^5 + 2*y^2 + 2).")
    print(f"This gives a natural factorization into h1(y) = y and h2(y) = y^5 + 2*y^2 + 2.")
    print(f"The resultant of h1(y) and h2(y) is h2(0).")
    print(f"h2(0) = 0**5 + 2*0**2 + 2 = {resultant_value}")
    print(f"The thickness is the 2-adic valuation of the resultant.")
    print(f"t = v({resultant_value})")
    print(f"Given v(2)=1, the thickness is {thickness}.")
    
solve()