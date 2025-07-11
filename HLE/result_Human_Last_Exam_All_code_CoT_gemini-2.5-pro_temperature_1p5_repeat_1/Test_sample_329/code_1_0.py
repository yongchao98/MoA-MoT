import numpy as np
from fractions import Fraction

def check_matrices():
    """
    Analyzes which of the given matrices are in the set P.
    P is the convex hull of matrices v*v^T for v in Z^2 \\ {(0,0)}.
    """
    
    # Define the matrices
    A = np.array([[0, 0], [0, 0]])
    B = np.array([[6, 4], [3, 7]])
    C = np.array([[1, -1/2], [-1/2, 1]])
    pi = np.pi
    D = np.array([[pi, 1], [1, pi**2]])
    E = np.array([[1, pi], [pi, 1]])
    F = np.array([[42, 0], [0, 0]])
    
    matrices = {'A': A, 'B': B, 'C': C, 'D': D, 'E': E, 'F': F}
    matrices_in_P = []

    print("Analyzing the properties of matrices in P = ConvexHull(v*v^T).")
    print("A matrix M in P must be symmetric, positive semidefinite, and have trace(M) >= 1.\n")

    for name, M in matrices.items():
        print(f"--- Checking Matrix {name} ---")
        print(M)
        
        # 1. Check for symmetry
        if not np.allclose(M, M.T):
            print(f"Reason: Matrix {name} is not symmetric ({M[0,1]} != {M[1,0]}).")
            print(f"Conclusion: {name} is not in P.\n")
            continue

        # 2. Check for positive semidefiniteness
        det = np.linalg.det(M)
        trace = np.trace(M)
        if det < -1e-9: # Allowing for small floating point errors
            print(f"Reason: Matrix {name} is not positive semidefinite (determinant = {det:.4f} < 0).")
            print(f"Conclusion: {name} is not in P.\n")
            continue
            
        # As it is symmetric, non-negative determinant and non-negative trace imply PSD for 2x2.
        # But we must check trace separately for the stronger condition tr(M) >= 1.
        if trace < -1e-9:
            print(f"Reason: Matrix {name} is not positive semidefinite (trace = {trace:.4f} < 0).")
            print(f"Conclusion: {name} is not in P.\n")
            continue

        # 3. Check trace >= 1 condition
        if trace < 1 - 1e-9:
            print(f"Reason: The trace of Matrix {name} is {trace:.4f}, which is less than 1.")
            print("The trace of any matrix in P must be >= 1.")
            print(f"Conclusion: {name} is not in P.\n")
            continue
        
        print(f"Matrix {name} is symmetric, positive semidefinite, and has trace >= 1.")

        # 4. Deeper analysis for matrices that pass the initial checks
        if name == 'C':
            print("Attempting to express C as a convex combination of matrices from S.")
            # Let's use v1 = [1, 1] and v2 = [1, -1]
            v1 = np.array([[1], [1]])
            S1 = v1 @ v1.T
            v2 = np.array([[1], [-1]])
            S2 = v2 @ v2.T
            # We want to find l1, l2 such that C = l1*S1 + l2*S2, with l1+l2=1, l1,l2>=0
            # From off-diagonal: -1/2 = l1*1 + l2*(-1) = l1 - l2
            # l1 + l2 = 1
            # Adding them: 2*l1 = 1/2 => l1 = 1/4. Then l2 = 3/4.
            l1 = Fraction(1, 4)
            l2 = Fraction(3, 4)
            print(f"C can be expressed as a convex combination: C = {l1}*S1 + {l2}*S2, where")
            print(f"S1 (from v=[1, 1]) is:\n{S1}")
            print(f"S2 (from v=[1, -1]) is:\n{S2}")
            print("The equation is:")
            print(f"[[1.  -0.5]\n [-0.5  1. ]] = {float(l1)} * [[{S1[0,0]}, {S1[0,1]}], [{S1[1,0]}, {S1[1,1]}]] + {float(l2)} * [[{S2[0,0]}, {S2[0,1]}], [{S2[1,0]}, {S2[1,1]}]]")
            print(f"Conclusion: C is in P.\n")
            matrices_in_P.append(name)
        elif name == 'D':
            print("Analysis for D:")
            print("If D were in P, it would have to lie in the affine space spanned by some vertices of P.")
            print("Since the vertices of P (matrices v*v^T) have integer entries, this affine space would be defined by a linear equation with rational coefficients: a*p + b*r + c*q + d = 0, where M = [[p, q], [q, r]].")
            print(f"For D, this means a*({pi:.4f}) + b*({pi**2:.4f}) + c*1 + d = 0 for some integers a,b,c,d (not all zero).")
            print("This would imply that pi is a root of a quadratic polynomial with integer coefficients.")
            print("However, pi is a transcendental number, meaning it is not a root of any non-zero polynomial with rational coefficients.")
            print("This is a contradiction.")
            print("Conclusion: D is not in P.\n")
        elif name == 'F':
            print("Attempting to express F as a convex combination of matrices from S.")
            # F = [[42, 0], [0, 0]]. For M=[[p,q],[q,r]] in P, p = sum(li*ai^2), r = sum(li*bi^2)
            # r = sum(li*bi^2) = 0 implies all bi=0 for li>0.
            # So F must be a convex combination of matrices from v=(a,0) form.
            # S_a = [[a^2, 0], [0, 0]].
            # We need 42 = sum(li*ai^2). This is possible since 42 is in ConvHull({1,4,9,16,25,36,49,...}) = [1, inf).
            # We can use a=6 and a=7. S_6 = [[36,0],[0,0]], S_7 = [[49,0],[0,0]].
            # 42 = l*36 + (1-l)*49 => 42 = 36l + 49 - 49l = 49 - 13l => 13l = 7 => l=7/13.
            l1 = Fraction(7, 13)
            l2 = Fraction(6, 13)
            v1 = np.array([[6],[0]])
            S1 = v1 @ v1.T
            v2 = np.array([[7],[0]])
            S2 = v2 @ v2.T
            print(f"F can be expressed as a convex combination: F = {l1}*S1 + {l2}*S2, where")
            print(f"S1 (from v=[6, 0]) is:\n{S1}")
            print(f"S2 (from v=[7, 0]) is:\n{S2}")
            print("The equation is:")
            print(f"[[42, 0], [0, 0]] = {l1.numerator}/{l1.denominator} * [[{S1[0,0]}, {S1[0,1]}], [{S1[1,0]}, {S1[1,1]}]] + {l2.numerator}/{l2.denominator} * [[{S2[0,0]}, {S2[0,1]}], [{S2[1,0]}, {S2[1,1]}]]")
            print(f"Conclusion: F is in P.\n")
            matrices_in_P.append(name)
            
    print("--- Summary ---")
    print(f"The matrices contained in P are: {matrices_in_P}")
    
    # Final answer formatting
    return matrices_in_P

if __name__ == '__main__':
    result = check_matrices()
    # The required format is a list of strings
    print(f"<<<[{','.join(sorted(result))}]>>>")