def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = 2  # Optimal number of disks to move to auxiliary peg for 4 disks
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Move remaining disks to target
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
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
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target)
    # Move nth disk to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Solve the problem with 4 disks, moving from Peg 2 to Peg 3
moves = hanoi_4_pegs(4, 2, 3, 1, 4)
for move in moves:
    print(move)