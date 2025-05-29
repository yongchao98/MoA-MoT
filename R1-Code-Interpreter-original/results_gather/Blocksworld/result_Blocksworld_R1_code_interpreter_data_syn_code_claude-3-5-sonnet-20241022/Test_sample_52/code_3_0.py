def print_solution():
    moves = [
        "Move C from 2 to 3",
        "Move E from 2 to 3",
        "Move G from 2 to 3",
        "Move B from 2 to 1",
        "Move F from 1 to 3",
        "Move J from 1 to 2",  # Temporarily move J
        "Move I from 1 to 3",  # Temporarily move I
        "Move H from 3 to 1",  # Move H to its position
        "Move D from 3 to 2",  # Move D to its final position
        "Move J from 2 to 1",  # Move J back
        "Move I from 3 to 1",  # Move I back
        "Move C from 3 to 1"   # Move C to its position
    ]
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()