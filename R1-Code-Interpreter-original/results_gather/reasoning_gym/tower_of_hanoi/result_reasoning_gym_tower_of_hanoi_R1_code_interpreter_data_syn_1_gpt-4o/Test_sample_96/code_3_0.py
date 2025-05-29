def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disks):
    if n == 0:
        return []
    if n == 1:
        return [(disks[0], source, target)]
    if n == 2:
        return [(disks[0], source, auxiliary1), (disks[1], source, target), (disks[0], auxiliary1, target)]
    else:
        moves = []
        # Move top n-2 disks to auxiliary1
        moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target, disks[:n-2])
        # Move next 2 disks to auxiliary2
        moves += hanoi_4_pegs(2, source, auxiliary2, target, auxiliary1, disks[n-2:])
        # Move 2 disks from auxiliary2 to target
        moves += hanoi_4_pegs(2, auxiliary2, target, source, auxiliary1, disks[n-2:])
        # Move n-2 disks from auxiliary1 to target
        moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2, disks[:n-2])
        return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(4, 3, 1, 2, 4, [1, 2, 3, 4])

# Format the moves
formatted_moves = [f"Move disk {disk} from Peg {move[1]} to Peg {move[2]}" for disk, move in moves]
print("\n".join(formatted_moves))