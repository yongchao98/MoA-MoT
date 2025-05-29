def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Use Frame-Stewart algorithm
    k = n - 1
    moves = []
    
    # Step 1: Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Step 2: Move the nth disk to target
    moves.append((n, source, target))
    
    # Step 3: Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Get the sequence of moves
moves = hanoi_4_pegs(3, 4, 1, 2, 3)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))