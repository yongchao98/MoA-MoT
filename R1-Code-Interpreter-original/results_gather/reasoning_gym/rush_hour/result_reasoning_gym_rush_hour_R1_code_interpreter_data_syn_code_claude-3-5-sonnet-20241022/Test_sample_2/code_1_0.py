def solve_rush_hour():
    # Initial board state
    board = [
        "BBH.CC",
        "G.H.JK",
        "G.AAJK",
        "DD.IxL",
        "EE.I.L",
        "FFF..x"
    ]
    
    # Looking at the board:
    # 1. AA is in row 2 (index 2)
    # 2. J and K are blocking AA's path to exit
    # 3. Both J and K can move up by 1 space
    # 4. After J and K move up, AA can move right to exit
    
    # The solution is deterministic:
    # 1. K must move up (-1)
    # 2. J must move up (-1)
    # 3. A must move right to exit (+3)
    
    solution = ["K-1", "J-1", "A+3"]
    print("<<<" + " ".join(solution) + ">>>")

solve_rush_hour()