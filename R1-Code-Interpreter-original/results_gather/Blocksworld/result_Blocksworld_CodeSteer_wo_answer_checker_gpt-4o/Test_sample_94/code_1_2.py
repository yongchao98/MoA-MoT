def blocksworld_solution():
    moves = [
        "Move G from 3 to 1",
        "Move H from 3 to 1",
        "Move C from 3 to 1",
        "Move E from 2 to 3",
        "Move I from 2 to 1",
        "Move D from 3 to 2"
    ]
    return moves

solution = blocksworld_solution()
print("<<<" + "\n".join(solution) + ">>>")