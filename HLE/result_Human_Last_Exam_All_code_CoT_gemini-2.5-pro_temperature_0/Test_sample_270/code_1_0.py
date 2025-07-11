import numpy as np

def solve_fdtc():
    """
    Calculates the fractional Dehn twist coefficient of (D_a * D_b)^9.
    """
    print("### Step 1: Represent Dehn twists as matrices in SL(2, Z)")
    print("The action of Dehn twists on the first homology of the torus can be represented by matrices.")
    print("Let 'a' and 'b' be curves generating the homology, with intersection number i(a, b) = 1.")
    # Matrix for a right-handed Dehn twist D_a about curve 'a'
    D_a = np.array([[1, 1], [0, 1]])
    # Matrix for a right-handed Dehn twist D_b about curve 'b'
    D_b = np.array([[1, 0], [-1, 1]])
    print(f"Matrix for D_a:\n{D_a}")
    print(f"Matrix for D_b:\n{D_b}\n")

    print("### Step 2: Compute the matrix for the composition M = D_a * D_b")
    M = D_a @ D_b
    print(f"Matrix M = D_a * D_b:\n{M}\n")

    print("### Step 3: Analyze the powers of M")
    print("We compute the matrices for (D_a * D_b)^3, (D_a * D_b)^6, and (D_a * D_b)^9.")
    M3 = np.linalg.matrix_power(M, 3)
    M6 = np.linalg.matrix_power(M, 6)
    M9 = np.linalg.matrix_power(M, 9)
    print(f"Matrix M^3:\n{M3}")
    print(f"Matrix M^6:\n{M6}")
    print(f"Matrix M^9:\n{M9}\n")

    print("### Step 4: Interpret the results in the mapping class group Mod(T^2_1)")
    print("The matrix M^6 is the identity matrix. This means the mapping class (D_a * D_b)^6 acts trivially on homology.")
    print("Such elements belong to the Torelli group and are known to be powers of the Dehn twist about the boundary curve, D_delta.")
    print("A standard relation in Mod(T^2_1) is (D_a * D_b)^6 = D_delta^(-1).")
    k = -1
    print(f"So, the power of the boundary twist is k = {k}.\n")

    print("### Step 5: Decompose the mapping class (D_a * D_b)^9")
    print("We write: (D_a * D_b)^9 = (D_a * D_b)^6 * (D_a * D_b)^3")
    print(f"Using the relation from Step 4, this becomes: D_delta^({k}) * (D_a * D_b)^3\n")

    print("### Step 6: Analyze the remaining term (D_a * D_b)^3")
    print("The matrix for (D_a * D_b)^3 is -I. This mapping class is the hyperelliptic involution.")
    print("The hyperelliptic involution is of finite order and can be represented by a homeomorphism that is the identity near the boundary.")
    print("The fractional Dehn twist coefficient (FDTC) of any such mapping class is 0.")
    c_rem = 0
    print(f"So, c((D_a * D_b)^3) = {c_rem}.\n")

    print("### Step 7: Calculate the final FDTC")
    print("The FDTC has the property that c(D_delta^k * psi) = k + c(psi).")
    final_coeff = k + c_rem
    print("For our element, the final equation is:")
    print(f"c((D_a * D_b)^9) = c(D_delta^({k})) + c((D_a * D_b)^3) = {k} + {c_rem} = {final_coeff}")
    
    print("\n--------------------------------------------------")
    print("The fractional Dehn twist coefficient is:")
    print(final_coeff)
    print("--------------------------------------------------")

solve_fdtc()