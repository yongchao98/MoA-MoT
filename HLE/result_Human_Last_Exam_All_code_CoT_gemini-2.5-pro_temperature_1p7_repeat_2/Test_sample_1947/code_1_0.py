def solve():
    """
    This function calculates the coefficients c_1, c_2, c_3, c_4, c_5 for the number of
    closed tree-like walks of length 6 in a simple graph X.

    The final expression is of the form:
    c_1*e + c_2*k + c_3*p + c_4*Sum(deg(v) choose 2) + c_5*Sum(deg(v) choose 3)

    Our derivation shows the coefficients to be:
    c_1 = 2
    c_2 = 0
    c_3 = 6
    c_4 = 12
    c_5 = 12
    """

    c1 = 2
    c2 = 0
    c3 = 6
    c4 = 12
    c5 = 12

    # Print the coefficients in the specified order
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")
    
    # Final answer format for the platform
    final_answer = f"<<<{c1}, {c2}, {c3}, {c4}, {c5}>>>"
    # This line is just for show in the final output block. The actual answer will be the last line.
    # print(final_answer)


# This function call will print the coefficients when executed
solve()
print("<<<2, 0, 6, 12, 12>>>")