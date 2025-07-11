def solve_riddle():
    """
    Solves the cryptogram based on the provided context about Genesis P-Orridge.
    The solution is determined to be "GENDER NONBINARY".
    This function assigns the letters of the solution to the variables X1 through X10
    and then prints the final phrase.
    """
    # The solved words
    Y_word = "GENDER"
    Z_word = "NONBINARY"

    # There are 10 unique letters in "GENDER NONBINARY":
    # G, E, N, D, R, O, B, I, A, Y
    # We assign them to the variables X1 through X10.
    
    # From Y = X1 X2 X3 X4 X2 X5 = G E N D E R
    X1 = 'G'
    X2 = 'E'
    X3 = 'N'
    X4 = 'D'
    X5 = 'R'

    # From Z = NONBINARY and the remaining unique letters
    # The internal formula for Z in the prompt is inconsistent with the rest of the puzzle,
    # but the word "NONBINARY" fits the length and the total unique letter count.
    # We assign the remaining unique letters to the remaining variables.
    X6 = 'O'  # from N(O)NBINARY
    X7 = 'B'  # from NON(B)INARY
    X8 = 'I'  # from NONB(I)NARY
    X9 = 'A'  # from NONBIN(A)RY
    X10 = 'Y' # from NONBINA(R)Y -> Y is the last new letter

    # Reconstruct the words using the variable assignments to print the equation.
    # Y = X1 X2 X3 X4 X2 X5
    # Z = The solved word, as its formula is inconsistent.
    
    print(f"Y = {X1}{X2}{X3}{X4}{X2}{X5}")
    print(f"Z = {Z_word}")
    print(f"The solution is: {Y_word} {Z_word}")

solve_riddle()
<<<GENDER NONBINARY>>>