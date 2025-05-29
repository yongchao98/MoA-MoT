def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(n, source, target)]
    
    # Find the optimal k
    k = n // 2
    moves = []
    
    # Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move remaining n-k disks to target
    moves.append((n, source, target))
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(6, 2, 3, 1, 4)

# Format the moves
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {tgt}" for disk, src, tgt in moves]

# Print the formatted moves
for move in formatted_moves:
    print(move)