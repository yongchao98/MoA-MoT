def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Choose k to minimize the number of moves
    k = 2  # For 5 disks, a common choice is to move 2 disks first
    moves = []
    
    # Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move remaining n-k disks to target
    moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Get the sequence of moves
moves = hanoi_4_pegs(5, 4, 1, 2, 3)

# Print the moves in the required format
for move in moves:
    print(f"Move disk {move[0]} from Peg {move[1]} to Peg {move[2]}")