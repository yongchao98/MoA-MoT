def print_moves():
    moves = [
        "Move F from 2 to 3",    # Clear A
        "Move A from 2 to 1",    # Place A in final position
        "Move F from 3 to 1",    # Place F in final position
        "Move H from 3 to 1",    # Place H in final position - Stack 1 complete
        "Move E from 3 to 2",    # Move E temporarily
        "Move I from 3 to 2",    # Move I temporarily
        "Move G from 2 to 1",    # Move G temporarily
        "Move E from 2 to 3",    # Move E to final stack
        "Move I from 2 to 3",    # Move I temporarily
        "Move C from 2 to 1",    # Move C temporarily
        "Move I from 3 to 2",    # Get I to final position
        "Move G from 1 to 2",    # Get G to final position
        "Move C from 1 to 2",    # Get C to final position - Stack 2 complete
        "Move D from 1 to 3"     # Complete Stack 3
    ]
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_moves()