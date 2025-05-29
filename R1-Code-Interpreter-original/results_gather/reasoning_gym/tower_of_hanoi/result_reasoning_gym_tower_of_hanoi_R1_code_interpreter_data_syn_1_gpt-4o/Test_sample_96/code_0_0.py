def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(source, target)]
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
formatted_moves = [f"Move disk {i+1} from Peg {move[0]} to Peg {move[1]}" for i, move in enumerate(moves)]
print("\n".join(formatted_moves))