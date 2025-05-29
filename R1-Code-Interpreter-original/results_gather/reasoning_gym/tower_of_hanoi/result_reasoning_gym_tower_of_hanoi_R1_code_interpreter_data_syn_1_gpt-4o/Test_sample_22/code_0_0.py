def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = n - 2  # Optimal k for 4 pegs and n disks
    # Step 1: Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Step 2: Move remaining disks to target
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    # Step 3: Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    # Move n-1 disks from source to auxiliary
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target)
    # Move nth disk from source to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Solve the problem for 4 disks, moving from Peg 3 to Peg 2
moves = hanoi_4_pegs(4, 3, 2, 1, 4)
for move in moves:
    print(move)