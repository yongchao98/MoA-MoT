def solve():
    """
    This function prints the computation of H_1(X, Z).
    The steps are:
    1. The space of nondegenerate lattices of unit area, X, is identified as the quotient
       X = (SL(2,R)/SL(2,Z)) / Z_2.
       Let Y = SL(2,R)/SL(2,Z). X = Y / Z_2.
    2. The first homology group of Y is H_1(Y, Z) = Z.
    3. The Z_2 action on Y induces an action on H_1(Y, Z) which corresponds to
       multiplication by -1.
    4. The homology of the quotient is the cokernel of the map (phi_* - id_*),
       where phi_* is the map on homology induced by the non-trivial element of Z_2.
       This map is (-1 - 1) = -2.
    5. So, H_1(X, Z) = coker(Z -> Z, z -> -2z) = Z / 2Z = Z_2.
    """
    
    # Representing the cyclic group of order 2.
    result = "Z_2"
    
    print(result)

solve()