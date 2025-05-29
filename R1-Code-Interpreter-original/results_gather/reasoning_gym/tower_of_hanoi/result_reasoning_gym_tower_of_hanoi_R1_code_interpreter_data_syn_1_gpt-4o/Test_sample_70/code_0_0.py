def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    # Move top n-2 disks to auxiliary1
    moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
    # Move the next two disks to target
    moves.append(f"Move disk {n-1} from Peg {source} to Peg {auxiliary2}")
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves.append(f"Move disk {n-1} from Peg {auxiliary2} to Peg {target}")
    # Move the n-2 disks from auxiliary1 to target
    moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves for 4 disks from Peg 2 to Peg 3
moves = hanoi_4_pegs(4, 2, 3, 1, 4)
print("\n".join(moves))