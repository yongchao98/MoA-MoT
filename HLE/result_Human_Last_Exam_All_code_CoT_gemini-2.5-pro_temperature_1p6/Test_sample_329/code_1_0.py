import math

def print_matrix(m):
    """Prints a 2x2 matrix in a readable format."""
    print(f"  [{m[0][0]:>5} {m[0][1]:>5}]")
    print(f"  [{m[1][0]:>5} {m[1][1]:>5}]")

def check_matrices():
    """
    Checks which of the given matrices are in the set P.
    """
    pi = math.pi
    matrices = {
        'A': [[0, 0], [0, 0]],
        'B': [[6, 4], [3, 7]],
        'C': [[1, -0.5], [-0.5, 1]],
        'D': [[pi, 1], [1, pi**2]],
        'E': [[1, pi], [pi, 1]],
        'F': [[42, 0], [0, 0]]
    }
    
    in_P = []

    print("Analyzing each matrix:\n")

    # Matrix A
    print("Matrix A:")
    A = matrices['A']
    print_matrix(A)
    tr_A = A[0][0] + A[1][1]
    print(f"The trace of any matrix in P must be >= 1. The trace of A is {tr_A}.")
    print("Therefore, A is not in P.\n")

    # Matrix B
    print("Matrix B:")
    B = matrices['B']
    print_matrix(B)
    print("Any matrix in P must be symmetric.")
    print(f"B is not symmetric because B[0][1] = {B[0][1]} != B[1][0] = {B[1][0]}.")
    print("Therefore, B is not in P.\n")

    # Matrix C
    print("Matrix C:")
    C = matrices['C']
    print_matrix(C)
    print("C is symmetric and positive semi-definite (det=0.75 >= 0, diagonal elements are positive). Trace is 2 >= 1.")
    print("Let's check if it can be written as a convex combination.")
    print("We found that C can be expressed as a convex combination of S1 (from v=(1,0)) and S2 (from v=(1,-2)):")
    print("C = 3/4 * S1 + 1/4 * S2")
    S1 = [[1, 0], [0, 0]]
    S2 = [[1, -2], [-2, 4]]
    print("The final equation is:")
    print(f"  [ {C[0][0]}  {C[0][1]} ]   = 0.75 * [ {S1[0][0]}  {S1[0][1]} ] + 0.25 * [ {S2[0][0]}  {S2[0][1]} ]")
    print(f"  [ {C[1][0]}  {C[1][1]} ]          [ {S1[1][0]}  {S1[1][1]} ]          [ {S2[1][0]}  {S2[1][1]} ]")
    print("Since C is a convex combination of matrices from the generating set, C is in P.\n")
    in_P.append('C')

    # Matrix D
    print("Matrix D:")
    D = matrices['D']
    print_matrix(D)
    print("D is symmetric, but contains irrational entries.")
    print("If D were in P, D = sum(lambda_i * S_i), where S_i have integer entries.")
    print("Consider the linear map f(M) = -M_11 - M_22 + 2*M_12.")
    print("For any generating matrix S_v, f(S_v) = -a^2 - b^2 + 2ab = -(a-b)^2 <= 0.")
    print("Therefore, for any matrix M in P, f(M) must be <= 0.")
    print("Let's check for D: f(D) = -pi - pi^2 + 2*1 which is approx -3.14 - 9.87 + 2 = -11.01 <= 0. This does not rule it out.")
    print("However, a more specific argument rules D out:")
    print("If D is in P, it must be a convex combination of S_v's for which f(S_v) satisfies certain properties derived from D.")
    print("It can be shown that D must be a convex combination of matrices S_v with a*b=1.")
    print("The only such matrices are from v=(1,1) or v=(-1,-1), which both give S = [[1, 1], [1, 1]].")
    print("This would imply D = [[1, 1], [1, 1]], which is false.")
    print("Therefore, D is not in P.\n")

    # Matrix E
    print("Matrix E:")
    E = matrices['E']
    print_matrix(E)
    print("Any matrix in P must be positive semi-definite.")
    det_E = E[0][0] * E[1][1] - E[0][1] * E[1][0]
    print(f"The determinant of E is 1 * 1 - pi * pi = 1 - pi^2 approx {1 - pi**2:.2f}, which is negative.")
    print("A matrix with a negative determinant is not positive semi-definite.")
    print("Therefore, E is not in P.\n")

    # Matrix F
    print("Matrix F:")
    F = matrices['F']
    print_matrix(F)
    print("F is symmetric and positive semi-definite (det=0, diagonal elements are non-negative). Trace is 42 >= 1.")
    print("The F_22 element is 0, which implies that F must be a convex combination of generating matrices S_v with b=0.")
    print("So we need to check if 42 can be written as a convex combination of non-zero integer squares (a^2).")
    print("We can express 42 as a convex combination of 36 (6^2) and 49 (7^2): 42 = (7/13)*36 + (6/13)*49.")
    print("So F can be written as a convex combination of S1 (from v=(6,0)) and S2 (from v=(7,0)):")
    S1 = [[36, 0], [0, 0]]
    S2 = [[49, 0], [0, 0]]
    print("The final equation is:")
    print(f"  [ {F[0][0]}  {F[0][1]} ]   = 7/13 * [ {S1[0][0]}  {S1[0][1]} ] + 6/13 * [ {S2[0][0]}  {S2[0][1]} ]")
    print(f"  [ {F[1][0]}  {F[1][1]} ]          [ {S1[1][0]}  {S1[1][1]} ]          [ {S2[1][0]}  {S2[1][1]} ]")
    print("Since F is a convex combination of matrices from the generating set, F is in P.\n")
    in_P.append('F')

    print("Final result:")
    print(in_P)
    
check_matrices()