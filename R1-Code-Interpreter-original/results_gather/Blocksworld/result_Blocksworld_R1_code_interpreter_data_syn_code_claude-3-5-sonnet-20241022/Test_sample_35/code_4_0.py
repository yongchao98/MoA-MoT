def print_solution():
    # Initial state: HGAC, E, F, BD
    # Goal state: ABEG, CF, DH
    
    moves = [
        "Move C from 1 to 2",    # HGAC, E, F, BD -> HGA, EC, F, BD
        "Move A from 1 to 3",    # HGA, EC, F, BD -> HG, EC, FA, BD
        "Move G from 1 to 4",    # HG, EC, FA, BD -> H, EC, FA, BDG
        "Move H from 1 to 3",    # H, EC, FA, BDG -> -, EC, FAH, BDG
        "Move D from 4 to 3",    # -, EC, FAH, BG -> -, EC, FAHD, BG
        "Move B from 4 to 1",    # -, EC, FAHD, G -> B, EC, FAHD, G
        "Move E from 2 to 1",    # B, C, FAHD, G -> BE, C, FAHD, G
        "Move G from 4 to 1",    # BE, C, FAHD, - -> BEG, C, FAHD
        "Move F from 3 to 2",    # BEG, C, AHD -> BEG, CF, AHD
        "Move A from 3 to 1",    # BEG, CF, HD -> ABEG, CF, HD
        "Move H from 3 to 3",    # ABEG, CF, D -> ABEG, CF, DH
    ]
    
    print('\n'.join(moves))

print_solution()