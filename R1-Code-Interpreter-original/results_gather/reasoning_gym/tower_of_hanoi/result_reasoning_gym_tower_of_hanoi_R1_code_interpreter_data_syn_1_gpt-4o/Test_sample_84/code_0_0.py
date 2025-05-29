def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(source, target)]
    
    # Find the optimal number of disks to move to an intermediate peg
    k = n - int((2 * n + 1) ** 0.5) + 1
    
    # Move k disks to auxiliary1
    moves = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move the remaining n-k disks to target
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(source, target)]
    
    moves = hanoi_3_pegs(n - 1, source, auxiliary, target)
    moves.append((source, target))
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(7, 1, 3, 2, 4)

# Format the output
formatted_moves = [f"Move disk {i+1} from Peg {move[0]} to Peg {move[1]}" for i, move in enumerate(moves)]
print("\n".join(formatted_moves))