def print_moves():
    moves = [
        "Move F from 2 to 3",    # Clear A
        "Move A from 2 to 1",    # Place A in final position
        "Move F from 3 to 1",    # Place F in final position
        "Move H from 3 to 1",    # Place H in final position
        "Move E from 3 to 2",    # Move E temporarily
        "Move I from 3 to 2",    # Place I in final position
        "Move E from 2 to 3",    # Move E to final position
        "Move D from 1 to 3",    # Move D to final position
        "Move C from 2 to 2",    # Place C in final position
        "Move G from 2 to 2"     # Place G in final position
    ]
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_moves()