def print_solution():
    moves = [
        "Move B from 2 to 1",     # G C B | D A E | I H F
        "Move E from 2 to 2",     # G C B | D A | I H F
        "Move A from 2 to 3",     # G C B | D | I H F A
        "Move D from 2 to 3",     # G C B | | I H F A D
        "Move C from 1 to 3",     # G B | | I H F A D C
        "Move F from 3 to 2",     # G B | F | I H A D C
        "Move G from 1 to 3",     # B | F | I H A D C G
        "Move F from 2 to 3",     # B | | I H A D C G F
        "Move I from 3 to 1",     # B I | | H A D C G F
        "Move H from 3 to 3",     # B I | E | A C D F G H
        "Move E from 2 to 2"      # B I | E | A C D F G H
    ]
    
    print('<<<')
    print('\n'.join(moves))
    print('>>>')

print_solution()