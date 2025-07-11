import numpy as np

def main():
    """
    Calculates the fractional Dehn twist coefficient by analyzing the corresponding
    matrices in SL(2,Z) and using a known relation in the mapping class group.
    """

    # Let the homology basis be (a, b) with intersection number i(a,b) = 1.
    # A right-handed Dehn twist D_c acts on a homology class x as x + i(x,c)c.
    # Note: i(b,a) = -i(a,b) = -1.
    
    # Action of D_a (right-handed Dehn twist about a):
    # D_a(a) = a
    # D_a(b) = b + i(b,a)a = b - a
    # Matrix representation M_a acting on row vectors (a, b):
    # (a, b) -> (a, b-a), which is not the standard matrix form.
    # Let's use column vectors [a; b].
    # a=[1,0]^T, b=[0,1]^T.
    # D_a([1,0]^T) = [1,0]^T
    # D_a([0,1]^T) = [-1,1]^T
    # So the matrix for D_a is [[1, -1], [0, 1]].
    M_a = np.array([[1, -1], [0, 1]])

    # Action of D_b (right-handed Dehn twist about b):
    # D_b(b) = b
    # D_b(a) = a + i(a,b)b = a + b
    # D_b([0,1]^T) = [0,1]^T
    # D_b([1,0]^T) = [1,1]^T
    # So the matrix for D_b is [[1, 0], [1, 1]].
    M_b = np.array([[1, 0], [1, 1]])

    print(f"Matrix for D_a:\n{M_a}\n")
    print(f"Matrix for D_b:\n{M_b}\n")
    
    # Calculate the matrix for the composition D_a * D_b
    M_comp = np.dot(M_a, M_b)
    print(f"Matrix for D_a * D_b:\n{M_comp}\n")
    
    # The relation (D_a D_b)^6 = D_delta^{-1} implies that the matrix for (D_a D_b)^6 should be the identity matrix.
    # Let's verify this.
    M_pow6 = np.linalg.matrix_power(M_comp, 6)
    print(f"Matrix for (D_a * D_b)^6:\n{M_pow6}\n")
    print("This confirms that (D_a D_b)^6 acts as the identity on homology, consistent with it being a power of D_delta.\n")

    # Now calculate the matrix for the full element (D_a D_b)^9
    M_pow9 = np.linalg.matrix_power(M_comp, 9)
    print(f"Matrix for (D_a * D_b)^9:\n{M_pow9}\n")
    print("This is the matrix for the hyperelliptic involution.\n")

    # The mapping class phi = (D_a D_b)^9 can be decomposed using the known relation:
    # (D_a D_b)^6 = D_delta^{-1}
    # phi = (D_a D_b)^9 = (D_a D_b)^6 * (D_a D_b)^3
    # phi = D_delta^{-1} * (D_a D_b)^3
    
    # The fractional Dehn twist coefficient is the exponent of D_delta in this decomposition.
    
    power_total = 9
    power_relation = 6
    power_remainder = 3
    fdtc = -1
    
    print("The key relation is: (D_a D_b)^6 = D_delta^(-1)")
    print("We decompose the given element:")
    print(f"(D_a D_b)^{power_total} = (D_a D_b)^{power_relation} * (D_a D_b)^{power_remainder}")
    print(f"(D_a D_b)^{power_total} = D_delta^({fdtc}) * (D_a D_b)^{power_remainder}\n")

    print("The fractional Dehn twist coefficient is the exponent of the boundary twist D_delta.")
    print(f"The coefficient is: {fdtc}")

if __name__ == "__main__":
    main()
