def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    # Find the optimal number of disks to move to the auxiliary peg
    k = n // 2
    
    # Move k disks to auxiliary1
    moves = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move the remaining n-k disks directly to the target
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = hanoi_3_pegs(n - 1, source, auxiliary, target)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Solve the problem for 7 disks, moving from Peg 3 to Peg 4
moves = hanoi_4_pegs(7, 3, 4, 1, 2)
print("\n".join(moves))