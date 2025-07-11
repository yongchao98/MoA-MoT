def solve_puzzle():
    """
    This puzzle follows a Fibonacci sequence pattern.
    The encoded words correspond to numbers in this sequence.
    """
    # The decoded numbers for the words are part of a Fibonacci sequence.
    # utpshtheas -> n1 = 8
    # fkuhu      -> n2 = 13
    # numhkatnatu-> n3 = 21
    n1 = 8
    n2 = 13
    n3 = n1 + n2

    # The final word, 'anhaeketn', decodes to the next number in the sequence.
    n4 = n2 + n3

    # The final equation is n2 + n3 = n4
    print(f"{n2} + {n3} = {n4}")

solve_puzzle()