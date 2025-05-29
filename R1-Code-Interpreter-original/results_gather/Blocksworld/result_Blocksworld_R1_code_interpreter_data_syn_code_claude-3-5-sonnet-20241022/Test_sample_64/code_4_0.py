def print_solution():
    moves = [
        "Move B from 2 to 1",     # G C B | D A E | I H F
        "Move E from 2 to 2",     # G C B | D A | I H F
        "Move A from 2 to 3",     # G C B | D | I H F A
        "Move D from 2 to 3",     # G C B | | I H F A D
        "Move C from 1 to 3",     # G B | | I H F A D C
        "Move B from 1 to 1",     # G B | | I H F A D C
        "Move I from 3 to 1",     # G B I | | H F A D C
        "Move H from 3 to 3",     # G B I | | F A D C H
        "Move G from 1 to 3",     # B I | | F A D C H G
        "Move E from 2 to 2"      # B I | E | A C D F G H
    ]
    
    print('<<<')
    print('\n'.join(moves))
    print('>>>')

print_solution()