def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    if n == 2:
        # Directly solve for 2 disks using 4 pegs
        moves.append(f"Move disk 1 from Peg {source} to Peg {auxiliary1}")
        moves.append(f"Move disk 2 from Peg {source} to Peg {target}")
        moves.append(f"Move disk 1 from Peg {auxiliary1} to Peg {target}")
    else:
        k = n - 2  # Optimal choice for 4 disks and 4 pegs

        # Step 1: Move top k disks to auxiliary1 using all 4 pegs
        moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
        
        # Step 2: Move remaining n-k disks to target using 3 pegs
        moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
        
        # Step 3: Move k disks from auxiliary1 to target using all 4 pegs
        moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    # Move n-1 disks from source to auxiliary
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    # Move nth disk from source to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves for 4 disks and 4 pegs
moves = hanoi_4_pegs(4, 4, 1, 2, 3)
print("\n".join(moves))