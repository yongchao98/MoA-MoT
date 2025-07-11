import math

def solve():
    """
    This function explains the reasoning to find the largest possible value of c.
    
    Let N be the number of planes in R^10.
    A point p is special if the vectors on all given planes through it span the whole R^10.
    The number of special points is always O(N^c). We want to find the largest possible value of c.

    Step 1: Analyze the condition for a point to be special.
    A plane in R^10 is a 2-dimensional affine subspace. Its associated vector space is a 2-dim linear subspace.
    For the vectors to span R^10, which is 10-dimensional, we need vectors from at least ceil(10 / 2) = 5 planes.
    So, a special point must be at the intersection of at least 5 planes.

    Step 2: Upper bound for c.
    Let's find an upper bound on the number of special points.
    A special point p must lie on at least 5 planes.
    The problem implies that the number of special points is finite for any configuration.
    As reasoned in the explanation, this forces the intersection of planes containing a special point to be that point itself.
    Let's count the number of intersections of 5 planes.
    The number of ways to choose 5 planes from N is C(N, 5).
    In a generic configuration designed to maximize intersections, any 5 planes intersect at a single point.
    Each such intersection point lies on at least 5 planes whose vector spaces (in general position) will span R^10. So these points are special.
    The number of such points is at most the number of ways to choose 5 planes from N, which is C(N, 5).
    C(N, 5) = N * (N-1) * (N-2) * (N-3) * (N-4) / (5 * 4 * 3 * 2 * 1)
    This is a polynomial in N of degree 5. So, the number of special points is O(N^5).
    This implies c <= 5.

    Step 3: Lower bound for c.
    We can construct a configuration of N planes that generates Omega(N^5) special points.
    Consider 5 sets of planes, each set containing N/5 parallel planes.
    Let the 5 vector spaces defining the orientation of each set of parallel planes be V1, V2, V3, V4, V5, chosen generically so that they span R^10.
    By arranging the planes appropriately (like a grid), the intersection of any 5 planes, one from each set, forms a unique point.
    The number of such intersection points is (N/5) * (N/5) * (N/5) * (N/5) * (N/5) = (N/5)^5 = N^5 / 3125.
    Each of these points is special because the 5 vector spaces span R^10.
    So, there exists a configuration with Omega(N^5) special points.
    This implies c >= 5.

    Step 4: Conclusion.
    From the upper bound (c <= 5) and the lower bound (c >= 5), we conclude that the largest possible value of c is 5.
    
    The equation for the number of special points in the upper bound argument is based on the number of combinations of 5 planes out of N.
    The number of combinations is given by the binomial coefficient "N choose 5".
    Let's write out the formula for N choose 5 to show the polynomial dependence.
    """
    
    c = 5
    
    print("Let N be the number of planes. The analysis suggests the maximum number of special points is bounded by the number of ways to choose 5 planes from N, which determines the order of the polynomial.")
    print("The number of combinations of 5 planes from N is given by the binomial coefficient C(N, 5).")
    print("C(N, 5) = N * (N-1) * (N-2) * (N-3) * (N-4) / (5 * 4 * 3 * 2 * 1)")
    print("This is a polynomial of degree 5 in N. For large N, this is approximately N^5 / 120.")
    print("Therefore, the number of special points is O(N^5).")
    print("A specific construction also shows that a number of special points of Omega(N^5) is achievable.")
    print(f"Combining the upper and lower bounds, the largest possible value of c is 5.")
    
    # We are asked to return the value of c.
    # No calculation is actually needed, the result is derived from the theoretical argument.
    
    print("\n---")
    print("Final Answer Derivation:")
    N = 'N'
    numerator = f"{N} * ({N}-1) * ({N}-2) * ({N}-3) * ({N}-4)"
    denominator = 5*4*3*2*1
    print(f"Number of special points is bounded by C(N, 5) = ({numerator}) / {denominator}")
    print("The dominant term is N^5, so the number of points is O(N^c) where c=5.")
    final_c = 5
    print(f"The largest possible value of c is {final_c}")

solve()