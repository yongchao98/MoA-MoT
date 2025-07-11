def solve_maze():
    """
    This function determines and prints the least dangerous path to the gold.
    The path is determined by analyzing the map provided by the user.
    """

    # The path is broken down into segments of movement in a single cardinal direction.
    # D = Down, U = Up, L = Left, R = Right

    # 1. Start at '@'. The only way out is down through the '+' door.
    move1 = "D"

    # 2. Inside the dragon's room. Move left to the far wall to maximize distance from the dragon.
    move2 = "L"

    # 3. Move down along the far wall.
    move3 = "D"
    
    # 4. Move left to exit the dragon's room through the far door and continue left down the hallway.
    move4 = "L"

    # 5. At the hallway junction, move up.
    move5 = "U"

    # 6. Move left through the final hallway, through the door to the gold's room, and to the gold 'g'.
    move6 = "L"

    # The final path is the sequence of these moves.
    final_path = move1 + move2 + move3 + move4 + move5 + move6
    
    print(final_path)

solve_maze()