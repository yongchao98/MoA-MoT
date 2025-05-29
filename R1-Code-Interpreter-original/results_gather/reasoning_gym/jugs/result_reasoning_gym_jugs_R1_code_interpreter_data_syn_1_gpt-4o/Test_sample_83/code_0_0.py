from collections import deque

def is_target_state(state, target):
    return target in state

def get_next_states(state):
    A, B, C = state
    next_states = []
    # Fill operations
    next_states.append((13, B, C))  # Fill A
    next_states.append((A, 13, C))  # Fill B
    next_states.append((A, B, 7))   # Fill C
    # Empty operations
    next_states.append((0, B, C))   # Empty A
    next_states.append((A, 0, C))   # Empty B
    next_states.append((A, B, 0))   # Empty C
    # Pour operations
    # Pour A -> B
    pour = min(A, 13 - B)
    next_states.append((A - pour, B + pour, C))
    # Pour A -> C
    pour = min(A, 7 - C)
    next_states.append((A - pour, B, C + pour))
    # Pour B -> A
    pour = min(B, 13 - A)
    next_states.append((A + pour, B - pour, C))
    # Pour B -> C
    pour = min(B, 7 - C)
    next_states.append((A, B - pour, C + pour))
    # Pour C -> A
    pour = min(C, 13 - A)
    next_states.append((A + pour, B, C - pour))
    # Pour C -> B
    pour = min(C, 13 - B)
    next_states.append((A, B + pour, C - pour))
    return next_states

def bfs_solve(target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_state(current_state, target):
            return path
        
        for next_state in get_next_states(current_state):
            if next_state not in visited:
                visited.add(next_state)
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

def format_moves(moves):
    formatted_moves = []
    for (start, end) in moves:
        if start[0] != end[0]:
            if end[0] == 13:
                formatted_moves.append("fill A")
            elif end[0] == 0:
                formatted_moves.append("empty A")
            else:
                formatted_moves.append("pour B->A" if start[1] != end[1] else "pour C->A")
        elif start[1] != end[1]:
            if end[1] == 13:
                formatted_moves.append("fill B")
            elif end[1] == 0:
                formatted_moves.append("empty B")
            else:
                formatted_moves.append("pour A->B" if start[0] != end[0] else "pour C->B")
        elif start[2] != end[2]:
            if end[2] == 7:
                formatted_moves.append("fill C")
            elif end[2] == 0:
                formatted_moves.append("empty C")
            else:
                formatted_moves.append("pour A->C" if start[0] != end[0] else "pour B->C")
    return formatted_moves

target = 11
moves = bfs_solve(target)
formatted_moves = format_moves(moves)
print(formatted_moves)