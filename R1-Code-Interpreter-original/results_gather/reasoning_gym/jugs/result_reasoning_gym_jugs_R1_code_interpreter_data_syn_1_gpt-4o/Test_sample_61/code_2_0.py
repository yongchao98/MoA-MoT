from collections import deque

def is_goal(state):
    return 7 in state

def get_next_states(state):
    a, b, c = state
    next_states = []
    # Fill operations
    next_states.append((3, b, c))  # Fill A
    next_states.append((a, 3, c))  # Fill B
    next_states.append((a, b, 11)) # Fill C
    # Empty operations
    next_states.append((0, b, c))  # Empty A
    next_states.append((a, 0, c))  # Empty B
    next_states.append((a, b, 0))  # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(a, 3 - b)
    next_states.append((a - pour, b + pour, c))
    # Pour A -> C
    pour = min(a, 11 - c)
    next_states.append((a - pour, b, c + pour))
    # Pour B -> A
    pour = min(b, 3 - a)
    next_states.append((a + pour, b - pour, c))
    # Pour B -> C
    pour = min(b, 11 - c)
    next_states.append((a, b - pour, c + pour))
    # Pour C -> A
    pour = min(c, 3 - a)
    next_states.append((a + pour, b, c - pour))
    # Pour C -> B
    pour = min(c, 3 - b)
    next_states.append((a, b + pour, c - pour))
    return next_states

def bfs():
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
        
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

def format_moves(moves):
    formatted_moves = []
    for (a1, b1, c1), (a2, b2, c2) in moves:
        if a1 != a2:
            if a2 > a1:
                formatted_moves.append("fill A")
            else:
                formatted_moves.append("empty A")
        elif b1 != b2:
            if b2 > b1:
                formatted_moves.append("fill B")
            else:
                formatted_moves.append("empty B")
        elif c1 != c2:
            if c2 > c1:
                formatted_moves.append("fill C")
            else:
                formatted_moves.append("empty C")
        else:
            if a1 != a2:
                formatted_moves.append("pour A->B")
            elif b1 != b2:
                formatted_moves.append("pour B->A")
            elif c1 != c2:
                formatted_moves.append("pour C->A")
    return formatted_moves

moves = bfs()
formatted_moves = format_moves(moves)
print(formatted_moves)