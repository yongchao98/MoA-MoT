import math

def solve():
    """
    Calculates the limit lim_{n->inf} [ln(1/p_n) / ln(n)].

    The problem asks for the escape probability p_n of a simple random walk from a discrete
    cube C_n = [0, 2n]^3, starting at S_0 = (n, 0, 0).

    1.  Interpretation: This is a problem about harmonic measure. The starting point
        S_0 = (n, 0, 0) is on an edge of the cube, where the faces y=0 and z=0 meet.
        "Escaping" is interpreted as the random walk hitting the other four "far" faces
        (x=0, x=2n, y=2n, z=2n) before hitting the two "near" faces (y=0, z=0).

    2.  Modeling: The probability p_n can be approximated by the value of a harmonic
        function u(x,y,z) at a point just inside the cube, e.g., (n, 1, 1). This
        function is the solution to the discrete Laplace equation with boundary conditions
        u=0 on the near faces and u=1 on the far faces.

    3.  Scaling Behavior: The continuum analogue is the Laplace equation for the potential
        V(x,y,z) in the cube. Near an edge where two grounded planes meet (a corner
        of co-dimension 2), the potential is proportional to the product of the
        distances to each plane. For a point like (n, 1, 1) in the cube [0, 2n]^3, the
        scaled distances to the faces y=0 and z=0 are proportional to 1/n.
        Thus, the escape probability p_n behaves as:
        p_n ~ C * (1/n) * (1/n) = C / n^2 for some constant C.

    4.  Limit Calculation: With p_n = C * n^(-2), we have:
        1/p_n = (1/C) * n^2
        ln(1/p_n) = ln(n^2 / C) = 2*ln(n) - ln(C)
        The limit is:
        lim_{n->inf} [ln(1/p_n) / ln(n)] = lim_{n->inf} [(2*ln(n) - ln(C)) / ln(n)]
                                       = lim_{n->inf} [2 - ln(C)/ln(n)] = 2.

    The final equation is based on the asymptotic behavior of p_n.
    p_n is approximately C / n^2.
    We are asked for lim_{n->inf} [ln(1/p_n) / ln(n)].
    Let's substitute the formula for p_n into the expression.
    Let alpha be the result.
    alpha = lim_{n->inf} [ln(n^2 / C) / ln(n)]
          = lim_{n->inf} [(2 * ln(n) - ln(C)) / ln(n)]
          = 2
    """
    
    alpha = 2
    
    # We are looking for the exponent in the power-law decay of p_n.
    # The derivation shows p_n is proportional to n^(-2).
    # The limit is the value of this exponent.
    
    print("The escape probability p_n scales as n^(-alpha).")
    print("p_n ~ C / n^2")
    print("Therefore, 1/p_n ~ (1/C) * n^2")
    print("ln(1/p_n) ~ ln((1/C) * n^2) = 2*ln(n) - ln(C)")
    print("The limit is lim_{n->inf} [ (2*ln(n) - ln(C)) / ln(n) ]")
    print("This simplifies to lim_{n->inf} [ 2 - ln(C)/ln(n) ] = 2.")
    print(f"The final answer is {alpha}")

solve()