def print_solution():
    # Initial state
    s1 = ['H', 'D']
    s2 = ['E', 'J', 'A', 'C', 'F']
    s3 = ['B', 'I', 'G', 'K']
    moves = []

    # Move F, C, A from stack 2 to clear top blocks
    moves.append("Move F from 2 to 3")
    moves.append("Move C from 2 to 3")
    moves.append("Move A from 2 to 3")
    moves.append("Move J from 2 to 3")
    moves.append("Move E from 2 to 1")

    # Build stack 1 (B E I J K)
    moves.append("Move K from 3 to 1")
    moves.append("Move I from 3 to 1")
    moves.append("Move B from 3 to 1")

    # Build stack 2 (A C D F)
    moves.append("Move D from 1 to 2")
    moves.append("Move F from 3 to 2")
    moves.append("Move C from 3 to 2")
    moves.append("Move A from 3 to 2")

    # Build stack 3 (G H)
    moves.append("Move G from 3 to 3")
    moves.append("Move H from 1 to 3")

    print("<<<")
    print('\n'.join(moves))
    print(">>>")

print_solution()