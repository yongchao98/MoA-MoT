def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Move top 1 disk to auxiliary1
    moves = hanoi_4_pegs(1, source, auxiliary1, auxiliary2, target)
    # Move remaining 2 disks to target using auxiliary2
    moves += hanoi_3_pegs(n-1, source, target, auxiliary2)
    # Move 1 disk from auxiliary1 to target
    moves += hanoi_4_pegs(1, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(3, 4, 1, 2, 3)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))