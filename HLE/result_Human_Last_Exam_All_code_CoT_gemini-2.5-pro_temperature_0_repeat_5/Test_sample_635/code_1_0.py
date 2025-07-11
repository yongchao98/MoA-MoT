import sympy

def solve_factorization():
    """
    This function calculates and displays the left coprime factorization
    of the given transfer function H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # The calculated left coprime factorization matrices D(s) and N(s)
    # H(s) = D(s)^-1 * N(s)
    D = sympy.Matrix([[s + 1, 0],
                      [1, s - 1]])

    N = sympy.Matrix([[s - 1, s + 1],
                      [1, 1]])

    # Output the final factorization in the required format
    print("A left coprime factorization of H(s) is given by H(s) = D(s)^-1 * N(s), where:")
    
    print("\nD(s) = ")
    print(f"[[{D[0,0]}, {D[0,1]}],")
    print(f" [{D[1,0]}, {D[1,1]}]]")

    print("\nN(s) = ")
    print(f"[[{N[0,0]}, {N[0,1]}],")
    print(f" [{N[1,0]}, {N[1,1]}]]")

    # --- Verification (optional) ---
    # This part verifies that the factorization is correct and coprime.
    # print("\n--- Verification ---")
    # H = sympy.Matrix([[(s - 1) / (s + 1), 1], [2 / (s**2 - 1), 0]])
    # D_inv_N = sympy.simplify(D.inv() * N)
    # assert sympy.simplify(D_inv_N - H) == sympy.zeros(2, 2)
    # print("Verification of H(s) = D(s)^-1 * N(s) successful.")
    
    # det_D = D.det()
    # roots = sympy.solve(det_D, s)
    # M = D.row_join(N)
    # is_coprime = True
    # for root in roots:
    #     if M.subs(s, root).rank() < D.rows:
    #         is_coprime = False
    #         break
    # assert is_coprime
    # print("Coprimeness verification successful.")

if __name__ == '__main__':
    solve_factorization()
