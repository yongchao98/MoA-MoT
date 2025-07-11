def solve_chess_puzzle():
    """
    This function prints the step-by-step move sequence for the correct answer.
    """
    # The chosen answer is C, as it's the only valid checkmating sequence.
    line_c = {
        1: ("Black", "1... Qg1+"),
        2: ("White", "2. Rxg1"),
        3: ("Black", "3. ... Nf2#")
    }

    print("The correct sequence leading to checkmate is:")
    print(f"{line_c[1][0]}'s move {line_c[1][1].split('...')[0]}: {line_c[1][1].split('... ')[1]}")
    print(f"{line_c[2][0]}'s move {line_c[2][1].split('.')[0]}: {line_c[2][1].split('. ')[1]}")
    print(f"{line_c[3][0]}'s move {line_c[3][1].split('...')[0]}: {line_c[3][1].split('... ')[1]}")

solve_chess_puzzle()