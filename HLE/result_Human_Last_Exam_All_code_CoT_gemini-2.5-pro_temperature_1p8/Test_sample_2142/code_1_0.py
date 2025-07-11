def solve_riddle():
    """
    Solves the P-Orridge riddle by deciphering the key and constructing the words.
    """
    # Step 1: Identify the keyword from the clue. The most relevant term for
    # "unification of masculine and feminine principles" is HERMAPHRODITE.
    keyword = "HERMAPHRODITE"

    # Step 2: Extract the 10 unique characters in order of appearance to form the key X.
    # The unique letters are H, E, R, M, A, P, O, D, I, T.
    X = ['H', 'E', 'R', 'M', 'A', 'P', 'O', 'D', 'I', 'T']
    
    # We map them to X1...X10 based on their index (X[0] is X1, etc.)
    X1, X2, X3, X4, X5, X6, X7, X8, X9, X10 = X

    # Step 3: Construct the word Y using the formula Y = X1X2X3X4X2X5.
    Y_chars = [X1, X2, X3, X4, X2, X5]
    Y = "".join(Y_chars)

    # Step 4: Construct the word Z using the formula Z = X3X6X1X7X8X9X10X2X6.
    Z_chars = [X3, X6, X1, X7, X8, X9, X10, X2, X6]
    Z = "".join(Z_chars)

    # Step 5: Print the construction of each word and the final answer.
    # The prompt asks to "output each number in the final equation".
    # This is interpreted as showing the component letters.
    
    print(f"The key characters X1-X10 are: {', '.join(X)}")
    print(f"Y = {' + '.join(Y_chars)} = {Y}")
    print(f"Z = {' + '.join(Z_chars)} = {Z}")
    
    final_answer = f"{Y} {Z}"
    print(f"\nThe commitment is to: {final_answer}")

solve_riddle()