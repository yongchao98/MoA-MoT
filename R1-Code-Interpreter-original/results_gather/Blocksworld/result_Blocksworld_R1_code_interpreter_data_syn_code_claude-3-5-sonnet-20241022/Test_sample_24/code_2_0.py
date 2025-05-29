def print_solution():
    # Initial state
    stack1 = ['I', 'D', 'B']
    stack2 = ['G', 'J', 'A', 'C', 'E']
    stack3 = ['H', 'F']
    
    # Pre-planned sequence of moves to reach the goal
    moves = [
        ('2', 'E', '1'),  # Move E to stack1
        ('3', 'F', '3'),  # Temporarily move F within stack3
        ('3', 'H', '1'),  # Move H to stack1
        ('2', 'C', '3'),  # Move C to stack3
        ('2', 'A', '2'),  # Keep A in stack2
        ('2', 'J', '3'),  # Move J to stack3
        ('1', 'B', '3'),  # Move B to stack3
        ('1', 'D', '3'),  # Move D to stack3
        ('1', 'I', '3'),  # Move I to stack3
        ('2', 'G', '2')   # Keep G in stack2
    ]
    
    print("<<<")
    for source, block, dest in moves:
        print(f"Move {block} from {source} to {dest}")
    print(">>>")

print_solution()