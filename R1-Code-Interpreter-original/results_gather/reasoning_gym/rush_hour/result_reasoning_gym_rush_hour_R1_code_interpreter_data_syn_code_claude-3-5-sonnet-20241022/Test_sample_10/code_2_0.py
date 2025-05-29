def print_solution():
    # Initial state:
    # .xBBCC
    # ..x.K.
    # G.AAK.
    # G.IJDD
    # H.IJ..
    # HEEFFF
    
    # Analysis:
    # 1. K needs to move down (+1) to clear space
    # 2. AA can then move right (+3) to reach exit
    
    solution = "K+1 A+3"
    print(f"<<<{solution}>>>")

print_solution()