import sympy

def v(n, p):
    """
    Computes the p-adic valuation of an integer n.
    Returns -1 for v(0) which we treat as infinity.
    """
    if n == 0:
        return -1 # Represents infinity
    count = 0
    while n % p == 0 and n != 0:
        count += 1
        n //= p
    return count

def main():
    """
    Calculates the thickness of the double point of the stable reduction.
    """
    # Step 1 & 2: Define the polynomial P(x)
    x = sympy.Symbol('x')
    P_x = 2*x**5 + 2*x**3 + 1
    
    # Step 3 & 4: Perform the change of variables x = 1 + 2y
    y = sympy.Symbol('y')
    # The double point corresponds to roots near x=1.
    # We substitute x = 1 + 2y to zoom in on this point.
    Q_y = P_x.subs(x, 1 + 2*y)
    
    # Step 5: Expand the transformed polynomial Q(y) and get its coefficients
    Q_y_expanded = sympy.expand(Q_y)
    
    # Extract coefficients of Q(y)
    # Q(y) = c5*y^5 + c4*y^4 + c3*y^3 + c2*y^2 + c1*y + c0
    coeffs = [Q_y_expanded.coeff(y, i) for i in range(6)]
    c0, c1, c2, c3, c4, c5 = [int(c) for c in coeffs]

    # Step 6: Analyze the Newton Polygon of Q(y) and calculate the thickness
    # We need the valuations of the coefficients of Q(y)
    v_coeffs = [v(c, 2) for c in coeffs]
    
    # The Newton Polygon for Q(y) has points (i, v(ci)).
    # (0, v(5)=0), (1, v(32)=5), (2, v(104)=3), ...
    # The lower convex hull contains a segment from (0, v(c0)) to (2, v(c2)).
    # This segment corresponds to two roots y1, y2.
    # The dominant terms for these roots are c2*y^2 + c0.
    # So y^2 is approximately -c0/c2.
    
    # The squared difference of the roots of Ay^2+B=0 is (2*sqrt(-B/A))^2 = -4B/A
    # A = c2, B = c0
    val_y_diff_sq = v(4, 2) + v(c0, 2) - v(c2, 2)
    
    # The thickness corresponds to the valuation of the squared difference of the x-roots.
    # e1 - e2 = (1+2*y1) - (1+2*y2) = 2*(y1 - y2)
    # (e1 - e2)^2 = 4 * (y1 - y2)^2
    thickness = v(4, 2) + val_y_diff_sq
    
    # Print out the results and the final equation
    print(f"The polynomial is P(x) = {P_x}")
    print(f"We perform the substitution x = 1 + 2y to get Q(y).")
    print(f"Q(y) = {Q_y_expanded}")
    print("\nThe coefficients of Q(y) are:")
    print(f"c0 = {c0}, c1 = {c1}, c2 = {c2}, c3 = {c3}, c4 = {c4}, c5 = {c5}")
    print("\nTheir 2-adic valuations are:")
    print(f"v(c0) = {v_coeffs[0]}, v(c1) = {v_coeffs[1]}, v(c2) = {v_coeffs[2]}, ...")
    print("\nThe Newton Polygon of Q(y) has a segment from (0,0) to (2,3).")
    print("This segment corresponds to two roots, y1 and y2.")
    print("The relevant terms are c2*y^2 and c0.")
    print(f"We approximate the equation for these roots as {c2}*y^2 + {c0} = 0")
    print("\nThe valuation of the squared difference of the y-roots is:")
    print(f"v((y1 - y2)^2) = v(4 * (-c0/c2)) = v(4) + v(c0) - v(c2)")
    print(f"v((y1 - y2)^2) = {v(4,2)} + {v(c0,2)} - {v(c2,2)} = {val_y_diff_sq}")

    print("\nThe thickness is the valuation of the squared difference of the original x-roots (e1, e2).")
    print("e1 - e2 = 2*(y1 - y2)")
    print("(e1 - e2)^2 = 4*(y1 - y2)^2")
    print("v((e1 - e2)^2) = v(4) + v((y1 - y2)^2)")
    print(f"Thickness = {v(4,2)} + {val_y_diff_sq} = {thickness}")

if __name__ == "__main__":
    main()
