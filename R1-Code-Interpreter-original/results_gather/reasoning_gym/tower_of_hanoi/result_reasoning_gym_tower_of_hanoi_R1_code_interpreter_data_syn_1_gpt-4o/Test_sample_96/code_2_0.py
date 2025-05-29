def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    if n == 2:
        return [(1, source, auxiliary1), (2, source, target), (1, auxiliary1, target)]
    else:
        moves = []
        # Move top n-2 disks to auxiliary1
        moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
        # Move next 2 disks to auxiliary2
        moves += hanoi_4_pegs(2, source, auxiliary2, target, auxiliary1)
        # Move 2 disks from auxiliary2 to target
        moves += hanoi_4_pegs(2, auxiliary2, target, source, auxiliary1)
        # Move n-2 disks from auxiliary1 to target
        moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
        return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(4, 3, 1, 2, 4)

# Format the moves
formatted_moves = [f"Move disk {disk} from Peg {move[1]} to Peg {move[2]}" for disk, move in enumerate(moves, start=1)]
print("\n".join(formatted_moves))