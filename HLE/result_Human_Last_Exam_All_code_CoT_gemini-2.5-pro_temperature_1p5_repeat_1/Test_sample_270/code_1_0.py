import math

def solve_fdtc():
    """
    This function calculates the fractional Dehn twist coefficient for the
    mapping class (Da o Db)^9 on a torus with one boundary component.
    """

    # The matrix representation of a right-handed Dehn twist about curve 'a'
    # which corresponds to the basis vector (1,0) in H_1.
    Da = [[1, 1], [0, 1]]

    # The matrix representation of a right-handed Dehn twist about curve 'b'
    # which corresponds to the basis vector (0,1) in H_1.
    Db = [[1, 0], [1, 1]]

    # A helper function for 2x2 matrix multiplication.
    # The composition f o g corresponds to the matrix product F * G.
    def mat_mul(A, B):
        return [
            [A[0][0]*B[0][0] + A[0][1]*B[1][0], A[0][0]*B[0][1] + A[0][1]*B[1][1]],
            [A[1][0]*B[0][0] + A[1][1]*B[1][0], A[1][0]*B[0][1] + A[1][1]*B[1][1]]
        ]

    # First, calculate the matrix for the single composition P = Da o Db.
    P = mat_mul(Da, Db)

    # Now, calculate the 9th power of P to find the matrix M for (Da o Db)^9.
    # We initialize M to P and multiply by P eight more times.
    M = P
    for _ in range(8):
        M = mat_mul(M, P)

    # Extract the elements p, q, r, s from the resulting matrix M = [[p, q], [r, s]]
    p = M[0][0]
    q = M[0][1]
    r = M[1][0]
    s = M[1][1]

    # The trace of M is p + s.
    trace_M = p + s

    # The formula for the FDTC is (p + s - 2*sign(p+s)) / r.
    # Since trace_M will be positive, sign(trace_M) is 1.
    numerator = trace_M - 2
    denominator = r
    
    # We now output the numbers in the final equation step-by-step.
    print(f"The matrix representation for (Da o Db)^9 is M = [[{p}, {q}], [{r}, {s}]].")
    print("The formula for the fractional Dehn twist coefficient (FDTC) is (p + s - 2*sign(p+s)) / r.")
    print("Plugging in the values from the matrix M:")
    print(f"FDTC = ({p} + {s} - 2*1) / {r}")
    print(f"FDTC = ({trace_M} - 2) / {denominator}")
    print(f"FDTC = {numerator} / {denominator}")

    # Simplify the resulting fraction by dividing by the greatest common divisor.
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print(f"The final simplified fraction is: {simplified_num}/{simplified_den}")


solve_fdtc()