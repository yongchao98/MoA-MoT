def print_moves(moves):
    result = ""
    for move in moves:
        result += f"Move {move[0]} from {move[1]} to {move[2]}\n"
    return result.strip()

def solve_blocksworld():
    moves = []
    
    # Clear space and get C D in stack1
    moves.append(('I', 4, 2))  # Move I to its final position in stack2
    moves.append(('H', 4, 1))  # Temporarily move H to stack1
    moves.append(('D', 4, 2))  # Temporarily move D to stack2
    moves.append(('G', 4, 1))  # Temporarily move G to stack1
    moves.append(('C', 4, 1))  # Move C to its final position in stack1
    moves.append(('G', 1, 3))  # Move G to stack3
    moves.append(('H', 1, 3))  # Move H to stack3
    moves.append(('D', 2, 1))  # Move D to its final position in stack1
    
    # Arrange stack2
    moves.append(('J', 3, 1))  # Temporarily move J
    moves.append(('F', 3, 2))  # Move F to its final position in stack2
    
    # Arrange stack3
    moves.append(('A', 3, 4))  # Temporarily move A
    moves.append(('B', 3, 1))  # Temporarily move B
    moves.append(('G', 3, 4))  # Temporarily move G
    moves.append(('H', 3, 4))  # Temporarily move H
    moves.append(('A', 4, 3))  # Move A to its final position
    moves.append(('B', 1, 3))  # Move B to its final position
    moves.append(('G', 4, 3))  # Move G to its final position
    moves.append(('H', 4, 3))  # Move H to its final position
    moves.append(('J', 1, 3))  # Move J to its final position
    
    print(print_moves(moves))

solve_blocksworld()