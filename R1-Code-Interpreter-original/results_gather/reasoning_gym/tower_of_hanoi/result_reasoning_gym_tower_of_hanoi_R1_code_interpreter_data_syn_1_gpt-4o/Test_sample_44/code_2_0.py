def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Find the optimal k
    k = n - 2
    moves = []
    
    # Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move remaining disks to target using auxiliary2
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Solve the problem for 5 disks and 4 pegs
moves = hanoi_4_pegs(5, 4, 1, 2, 3)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))