def tower_of_hanoi_4_pegs():
    moves = [
        "Move disk 1 from Peg 1 to Peg 3",
        "Move disk 2 from Peg 1 to Peg 4",
        "Move disk 1 from Peg 3 to Peg 2",
        "Move disk 3 from Peg 1 to Peg 3",
        "Move disk 1 from Peg 2 to Peg 4",
        "Move disk 2 from Peg 4 to Peg 2",
        "Move disk 1 from Peg 4 to Peg 2",
        "Move disk 3 from Peg 3 to Peg 2"
    ]
    for move in moves:
        print(move)

tower_of_hanoi_4_pegs()