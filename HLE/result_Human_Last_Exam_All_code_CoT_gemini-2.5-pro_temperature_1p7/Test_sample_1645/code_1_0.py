def find_smallest_n():
    """
    This function determines the smallest non-negative integer n such that the property (Rn)
    is not preserved by the completion of a noetherian local ring.

    The property (Rn) for a noetherian ring A means that for every prime ideal p
    with height(p) <= n, the localization A_p is a regular local ring.

    The analysis proceeds step-by-step for n = 0, 1, ...
    """

    # Case n = 0:
    # (R0) is equivalent to the ring A being reduced.
    # For a noetherian local ring A, A is reduced if and only if its completion Â is reduced.
    # Thus, the property (R0) is always preserved under completion.
    # So the answer is not 0.
    is_R0_preserved = True

    # Case n = 1:
    # For a noetherian domain, (R1) is equivalent to the ring being normal
    # (by Serre's criterion for normality).
    # A famous counterexample by Nagata shows there exists a normal noetherian local domain A
    # (which thus satisfies (R1)) for which its completion Â is not normal
    # (and thus does not satisfy (R1)).
    # Therefore, the property (R1) is not always preserved.
    is_R1_preserved = False
    
    # We are looking for the smallest n for which the property is not preserved.
    # Since it's preserved for n=0 and not for n=1, the smallest such integer is 1.
    smallest_n = -1
    if is_R0_preserved and not is_R1_preserved:
        smallest_n = 1

    # The problem asks for an integer n, so we will output it directly.
    # Based on the established theory of Nagata rings:
    n = 1
    
    # We don't have an equation, so we print the resulting number directly.
    print(n)

find_smallest_n()