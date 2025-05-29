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
        # Move top n-2 disks to auxiliary1 using all 4 pegs
        moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
        
        # Move the next two disks to the target using auxiliary2 as temporary storage
        moves.append(f"Move disk {n-1} from Peg {source} to Peg {auxiliary2}")
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        moves.append(f"Move disk {n-1} from Peg {auxiliary2} to Peg {target}")
        
        # Move the n-2 disks from auxiliary1 to target using all 4 pegs
        moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves for 4 disks and 4 pegs
moves = hanoi_4_pegs(4, 4, 1, 2, 3)
print("\n".join(moves))