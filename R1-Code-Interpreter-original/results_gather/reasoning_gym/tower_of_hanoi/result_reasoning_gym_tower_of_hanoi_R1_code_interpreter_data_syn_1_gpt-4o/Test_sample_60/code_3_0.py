def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Optimal k for 6 disks and 4 pegs is known to be 3
    k = 3
    
    # Step 1: Move top k disks to auxiliary1 using all pegs
    moves = hanoi_4_pegs(k, source, auxiliary1, target, auxiliary2)
    
    # Step 2: Move remaining n-k disks to target using 3 pegs
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    
    # Step 3: Move k disks from auxiliary1 to target using all pegs
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = hanoi_3_pegs(n - 1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(6, 1, 3, 2, 4)

# Format the moves as required
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Print the formatted moves
for move in formatted_moves:
    print(move)