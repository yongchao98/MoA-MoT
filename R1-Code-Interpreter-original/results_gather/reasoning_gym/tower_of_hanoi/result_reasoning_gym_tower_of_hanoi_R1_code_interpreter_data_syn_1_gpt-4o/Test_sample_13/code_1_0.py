def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = n - 2  # For 5 disks, a common choice is to move 3 disks first
    
    # Step 1: Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Step 2: Move remaining disks to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    
    # Step 3: Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Solve the problem for 5 disks, moving from Peg 4 to Peg 2
moves = hanoi_4_pegs(5, 4, 2, 1, 3)
print("\n".join(moves))