def print_solution():
    # Initial state: G C | D A E B | I H F
    # Goal state:   B I | E       | A C D F G H
    
    # Direct solution with minimal moves
    moves = [
        "Move B from 2 to 1",    # G C B | D A E | I H F
        "Move E from 2 to 2",    # G C B | D A | I H F
        "Move A from 2 to 3",    # G C B | D | I H F A
        "Move D from 2 to 3",    # G C B | | I H F A D
        "Move C from 1 to 3",    # G B | | I H F A D C
        "Move B from 1 to 1",    # G | | I H F A D C B
        "Move G from 1 to 3",    # | | I H F A D C B G
        "Move I from 3 to 1",    # I | H F A D C B G
        "Move H from 3 to 3"     # B I | E | A C D F G H
    ]
    
    print('<<<')
    print('\n'.join(moves))
    print('>>>')

print_solution()