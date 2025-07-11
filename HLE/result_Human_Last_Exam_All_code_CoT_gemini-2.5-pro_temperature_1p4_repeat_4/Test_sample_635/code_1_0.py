def solve_left_coprime_factorization():
    """
    This function prints the matrices for a left coprime factorization
    of the given transfer function H(s).
    The factorization is in the form H(s) = D(s)^-1 * N(s).
    """

    print("The left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("")

    # Display the D(s) matrix
    print("D(s) =")
    print("[[ s + 1,     0 ]]")
    print("[[     1, s - 1 ]]")
    print("")

    # Display the N(s) matrix
    print("N(s) =")
    print("[[ s - 1, s + 1 ]]")
    print("[[     1,     1 ]]")
    print("")
    
    # Display the final equation using the matrices
    print("Therefore, the final equation is:")
    # The string formatting below aligns the matrices for clear visualization.
    # The numbers in the string represent the characters for each matrix row.
    print(f"{'H(s) = [[ s + 1,     0 ]]'} ^-1  *  {'[[ s - 1, s + 1 ]]'}")
    print(f"{'       [[     1, s - 1 ]]'}          {'[[     1,     1 ]]'}")

solve_left_coprime_factorization()