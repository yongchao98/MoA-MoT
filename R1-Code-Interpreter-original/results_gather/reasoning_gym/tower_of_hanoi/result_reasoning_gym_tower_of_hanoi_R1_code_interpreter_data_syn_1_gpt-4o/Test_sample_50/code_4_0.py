def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    # Find the optimal k
    k = n - 1
    for i in range(1, n):
        if 2 * len(hanoi_3_pegs(i, source, auxiliary1, auxiliary2)) + len(hanoi_3_pegs(n-i, source, target, auxiliary2)) < 2 * len(hanoi_3_pegs(k, source, auxiliary1, auxiliary2)) + len(hanoi_3_pegs(n-k, source, target, auxiliary2)):
            k = i
    
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
    # Move n-1 disks to auxiliary
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    # Move nth disk to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Solve the problem for 4 disks, moving from Peg 3 to Peg 1
moves = hanoi_4_pegs(4, 3, 1, 2, 4)
print("\n".join(moves))