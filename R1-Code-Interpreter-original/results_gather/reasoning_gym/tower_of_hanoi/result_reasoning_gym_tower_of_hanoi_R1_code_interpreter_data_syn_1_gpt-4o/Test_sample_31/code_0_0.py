def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = n - 2  # Optimal k for 6 disks and 4 pegs is 4
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Move remaining n-k disks to target using auxiliary2
    moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Solve the problem for 6 disks, moving from Peg 3 to Peg 2
moves = hanoi_4_pegs(6, 3, 2, 1, 4)
print("\n".join(moves))