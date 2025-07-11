def solve_maze():
    """
    This function determines the least dangerous path from the adventurer '@' to the gold '$g'.

    The path is determined by the following logic:
    1.  The primary danger, a dragon 'D', must be avoided.
    2.  The adventurer starts in a room with only one exit, which is downwards.
    3.  This leads to a corridor with the dragon below. To avoid it, the path must go up.
    4.  From the top of the map, the path moves left.
    5.  Then, the path moves down to the correct level for the gold room's entrance.
    6.  Finally, the path moves left into the room to get the gold.
    """
    
    # Sequence of moves: Down, Up, Left, Down, Left
    path = "DULDL"
    
    print(path)

solve_maze()