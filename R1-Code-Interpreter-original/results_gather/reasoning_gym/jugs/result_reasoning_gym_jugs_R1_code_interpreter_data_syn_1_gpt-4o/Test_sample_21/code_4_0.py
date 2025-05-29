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
    if A > 0 and B < capB:
        pour_amount = min(A, capB - B)
        moves.append((A - pour_amount, B + pour_amount, C))
    # Pour A -> C
    if A > 0 and C < capC:
        pour_amount = min(A, capC - C)
        moves.append((A - pour_amount, B, C + pour_amount))
    # Pour B -> A
    if B > 0 and A < capA:
        pour_amount = min(B, capA - A)
        moves.append((A + pour_amount, B - pour_amount, C))
    # Pour B -> C
    if B > 0 and C < capC:
        pour_amount = min(B, capC - C)
        moves.append((A, B - pour_amount, C + pour_amount))
    # Pour C -> A
    if C > 0 and A < capA:
        pour_amount = min(C, capA - A)
        moves.append((A + pour_amount, B, C - pour_amount))
    # Pour C -> B
    if C > 0 and B < capB:
        pour_amount = min(C, capB - B)
        moves.append((A, B + pour_amount, C - pour_amount))
    
    return moves

def find_solution(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_reached(current_state, target):
            return path + [current_state]
        
        for next_state in get_possible_moves(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + [next_state]))
    
    return None

def convert_states_to_moves(states):
    moves = []
    for i in range(1, len(states)):
        prev_state = states[i - 1]
        curr_state = states[i]
        
        if prev_state[0] != curr_state[0]:
            if curr_state[0] == 4:
                moves.append("fill A")
            elif curr_state[0] == 0:
                moves.append("empty A")
            else:
                if prev_state[1] != curr_state[1]:
                    moves.append("pour B->A")
                else:
                    moves.append("pour C->A")
        
        elif prev_state[1] != curr_state[1]:
            if curr_state[1] == 9:
                moves.append("fill B")
            elif curr_state[1] == 0:
                moves.append("empty B")
            else:
                if prev_state[0] != curr_state[0]:
                    moves.append("pour A->B")
                else:
                    moves.append("pour C->B")
        
        elif prev_state[2] != curr_state[2]:
            if curr_state[2] == 9:
                moves.append("fill C")
            elif curr_state[2] == 0:
                moves.append("empty C")
            else:
                if prev_state[0] != curr_state[0]:
                    moves.append("pour A->C")
                else:
                    moves.append("pour B->C")
    
    return moves

capacities = (4, 9, 9)
target = 2
solution = find_solution(capacities, target)
moves = convert_states_to_moves(solution)
print(moves)