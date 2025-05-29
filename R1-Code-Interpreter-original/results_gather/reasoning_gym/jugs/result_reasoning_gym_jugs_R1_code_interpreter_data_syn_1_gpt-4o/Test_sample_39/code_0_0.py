from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    moves = []
    
    # Fill operations
    if A < capA:
        moves.append((capA, B, C))
    if B < capB:
        moves.append((A, capB, C))
    if C < capC:
        moves.append((A, B, capC))
    
    # Empty operations
    if A > 0:
        moves.append((0, B, C))
    if B > 0:
        moves.append((A, 0, C))
    if C > 0:
        moves.append((A, B, 0))
    
    # Pour operations
    # Pour A -> B
    pour = min(A, capB - B)
    moves.append((A - pour, B + pour, C))
    
    # Pour A -> C
    pour = min(A, capC - C)
    moves.append((A - pour, B, C + pour))
    
    # Pour B -> A
    pour = min(B, capA - A)
    moves.append((A + pour, B - pour, C))
    
    # Pour B -> C
    pour = min(B, capC - C)
    moves.append((A, B - pour, C + pour))
    
    # Pour C -> A
    pour = min(C, capA - A)
    moves.append((A + pour, B, C - pour))
    
    # Pour C -> B
    pour = min(C, capB - B)
    moves.append((A, B + pour, C - pour))
    
    return moves

def find_solution(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_reached(current_state, target):
            return path
        
        if current_state in visited:
            continue
        
        visited.add(current_state)
        
        for next_state in get_possible_moves(current_state, capacities):
            if next_state not in visited:
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

capacities = (13, 12, 6)
target = 3
solution = find_solution(capacities, target)
print(solution)