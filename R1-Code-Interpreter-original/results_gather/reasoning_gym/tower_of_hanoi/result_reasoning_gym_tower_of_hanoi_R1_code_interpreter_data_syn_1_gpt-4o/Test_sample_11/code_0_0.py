def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    # Move n-2 disks from source to auxiliary1 using auxiliary2 and target
    moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
    # Move the (n-1)th disk to auxiliary2
    moves.append(f"Move disk {n-1} from Peg {source} to Peg {auxiliary2}")
    # Move the nth disk to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move the (n-1)th disk from auxiliary2 to target
    moves.append(f"Move disk {n-1} from Peg {auxiliary2} to Peg {target}")
    # Move n-2 disks from auxiliary1 to target using source and auxiliary2
    moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves for 3 disks from Peg 2 to Peg 4
moves = hanoi_4_pegs(3, 2, 4, 1, 3)
print("\n".join(moves))