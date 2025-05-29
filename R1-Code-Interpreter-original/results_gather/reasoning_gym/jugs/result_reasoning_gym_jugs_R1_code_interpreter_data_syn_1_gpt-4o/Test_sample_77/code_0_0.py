from collections import deque

# Define the capacities of the jugs
capacities = (10, 3, 13)
target = 8

# Define the initial state (all jugs are empty)
initial_state = (0, 0, 0)

# Define a function to get the next states from a given state
def get_next_states(state):
    a, b, c = state
    next_states = []
    
    # Fill operations
    next_states.append((capacities[0], b, c))  # Fill A
    next_states.append((a, capacities[1], c))  # Fill B
    next_states.append((a, b, capacities[2]))  # Fill C
    
    # Empty operations
    next_states.append((0, b, c))  # Empty A
    next_states.append((a, 0, c))  # Empty B
    next_states.append((a, b, 0))  # Empty C
    
    # Pour operations
    # Pour A -> B
    pour = min(a, capacities[1] - b)
    next_states.append((a - pour, b + pour, c))
    
    # Pour A -> C
    pour = min(a, capacities[2] - c)
    next_states.append((a - pour, b, c + pour))
    
    # Pour B -> A
    pour = min(b, capacities[0] - a)
    next_states.append((a + pour, b - pour, c))
    
    # Pour B -> C
    pour = min(b, capacities[2] - c)
    next_states.append((a, b - pour, c + pour))
    
    # Pour C -> A
    pour = min(c, capacities[0] - a)
    next_states.append((a + pour, b, c - pour))
    
    # Pour C -> B
    pour = min(c, capacities[1] - b)
    next_states.append((a, b + pour, c - pour))
    
    return next_states

# Perform BFS to find the solution
def bfs():
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we have reached the target
        if target in current_state:
            return path
        
        # Get the next states
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                # Determine the move that led to this state
                move = determine_move(current_state, next_state)
                queue.append((next_state, path + [move]))

# Determine the move that led to the next state
def determine_move(current, next):
    a1, b1, c1 = current
    a2, b2, c2 = next
    
    if a1 != a2:
        if a2 == capacities[0]:
            return 'fill A'
        elif a2 == 0:
            return 'empty A'
        elif a1 > a2:
            return 'pour A->' + ('B' if b1 != b2 else 'C')
        else:
            return 'pour ' + ('B' if b1 != b2 else 'C') + '->A'
    
    if b1 != b2:
        if b2 == capacities[1]:
            return 'fill B'
        elif b2 == 0:
            return 'empty B'
        elif b1 > b2:
            return 'pour B->' + ('A' if a1 != a2 else 'C')
        else:
            return 'pour ' + ('A' if a1 != a2 else 'C') + '->B'
    
    if c1 != c2:
        if c2 == capacities[2]:
            return 'fill C'
        elif c2 == 0:
            return 'empty C'
        elif c1 > c2:
            return 'pour C->' + ('A' if a1 != a2 else 'B')
        else:
            return 'pour ' + ('A' if a1 != a2 else 'B') + '->C'

# Find the solution
solution = bfs()
print(solution)