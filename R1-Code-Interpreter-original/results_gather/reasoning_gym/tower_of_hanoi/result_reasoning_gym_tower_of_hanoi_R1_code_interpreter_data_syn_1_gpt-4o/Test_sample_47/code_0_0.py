def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    k = n - 2  # Optimal k for 4 pegs and n disks
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Move remaining n-k disks to target using 3 pegs
    moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    # Move n-1 disks to auxiliary
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    # Move nth disk to target
    moves.append((n, source, target))
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(4, 3, 4, 1, 2)

# Format the moves
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Print the formatted moves
print("\n".join(formatted_moves))