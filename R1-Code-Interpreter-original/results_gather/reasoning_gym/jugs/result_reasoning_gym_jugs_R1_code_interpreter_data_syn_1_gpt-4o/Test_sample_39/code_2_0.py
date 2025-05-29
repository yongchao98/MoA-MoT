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

def convert_to_moves(state_transitions):
    moves = []
    for (prev_state, next_state) in state_transitions:
        A1, B1, C1 = prev_state
        A2, B2, C2 = next_state
        
        if A1 != A2:
            if A2 == 13:
                moves.append("fill A")
            elif A2 == 0:
                moves.append("empty A")
            else:
                if A1 > A2:
                    moves.append("pour A->B" if B1 != B2 else "pour A->C")
                else:
                    moves.append("pour B->A" if B1 != B2 else "pour C->A")
        
        elif B1 != B2:
            if B2 == 12:
                moves.append("fill B")
            elif B2 == 0:
                moves.append("empty B")
            else:
                if B1 > B2:
                    moves.append("pour B->A" if A1 != A2 else "pour B->C")
                else:
                    moves.append("pour A->B" if A1 != A2 else "pour C->B")
        
        elif C1 != C2:
            if C2 == 6:
                moves.append("fill C")
            elif C2 == 0:
                moves.append("empty C")
            else:
                if C1 > C2:
                    moves.append("pour C->A" if A1 != A2 else "pour C->B")
                else:
                    moves.append("pour A->C" if A1 != A2 else "pour B->C")
    
    return moves

capacities = (13, 12, 6)
target = 3
state_transitions = find_solution(capacities, target)
moves = convert_to_moves(state_transitions)
print(moves)