def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(n, source, target)]
    
    # Find the optimal k
    k = n - 1
    moves = []
    
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move the remaining n-k disks to target
    moves += hanoi_4_pegs(n-k, source, target, auxiliary2, auxiliary1)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Get the sequence of moves for 5 disks
moves = hanoi_4_pegs(5, 1, 2, 3, 4)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))