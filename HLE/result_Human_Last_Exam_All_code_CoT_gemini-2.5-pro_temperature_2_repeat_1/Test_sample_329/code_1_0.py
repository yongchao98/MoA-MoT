def solve():
    """
    This function determines which of the given matrices are in the set P.

    The set P is the convex hull of matrices v * v^T where v is a non-zero 2D integer vector.
    A matrix M is in P if it's a convex combination of such generating matrices.

    Properties of matrices in P:
    1. Symmetric
    2. Positive semidefinite
    3. Has a positive trace

    Matrix-by-matrix analysis:
    A = [[0,0],[0,0]]: Fails property 3 (trace is 0).
    B = [[6,4],[3,7]]: Fails property 1 (not symmetric).
    C = [[1,-1/2],[-1/2,1]]: All basic properties hold. It can be written as a convex combination of S_v matrices.
       Let's check C = l1 * [[1,1],[1,1]] + l2 * [[1,-1],[-1,1]], where the generators are for v=(1,1) and v=(1,-1).
       l1 + l2 = 1 and l1 - l2 = -1/2. Solving gives l1=1/4, l2=3/4. Both are positive. So C is in P.
    D = [[pi,1],[1,pi^2]]: The relation D_22 = D_11^2 implies, if D is in P, that pi must be rational, which is false. So D is not in P.
    E = [[1,pi],[pi,1]]: Fails property 2 (determinant is 1 - pi^2 < 0, so not positive semidefinite).
    F = [[42,0],[0,0]]: All basic properties hold. It can be a convex combination of matrices for v=(a,0).
       We need to make 42 from integer squares. For instance, using 36 (from v=(6,0)) and 49 (from v=(7,0)).
       42 = l1 * 36 + l2 * 49, with l1+l2=1. This gives l1=7/13 and l2=6/13. So F is in P.

    Thus, the matrices in P are C and F.
    """
    
    result = ['C', 'F']
    print(result)

solve()