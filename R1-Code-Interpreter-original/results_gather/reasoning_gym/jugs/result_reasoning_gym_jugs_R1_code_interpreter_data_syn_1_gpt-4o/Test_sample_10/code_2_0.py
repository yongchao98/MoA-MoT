from collections import deque

def is_target_reached(state, target):
    return target in state

def get_possible_moves(state, capacities):
    A, B, C = state
    capA, capB, capC = capacities
    moves = []
    
    # Fill operations
    if A < capA:
        moves.append((capA, B, C))  # Fill A
    if B < capB:
        moves.append((A, capB, C))  # Fill B
    if C < capC:
        moves.append((A, B, capC))  # Fill C
    
    # Empty operations
    if A > 0:
        moves.append((0, B, C))  # Empty A
    if B > 0:
        moves.append((A, 0, C))  # Empty B
    if C > 0:
        moves.append((A, B, 0))  # Empty C
    
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
            return path
        
        for next_state in get_possible_moves(current_state, capacities):
            if next_state not in visited:
                visited.add(next_state)
                move = (current_state, next_state)
                queue.append((next_state, path + [move]))
    
    return None

capacities = (13, 13, 4)
target = 10
solution = find_solution(capacities, target)

# Convert the solution to a JSON-parsable list of moves
moves = []
for (start, end) in solution:
    if start[0] != end[0]:
        if end[0] == capacities[0]:
            moves.append("fill A")
        elif end[0] == 0:
            moves.append("empty A")
        else:
            if start[1] != end[1]:
                moves.append("pour B->A")
            else:
                moves.append("pour C->A")
    elif start[1] != end[1]:
        if end[1] == capacities[1]:
            moves.append("fill B")
        elif end[1] == 0:
            moves.append("empty B")
        else:
            if start[0] != end[0]:
                moves.append("pour A->B")
            else:
                moves.append("pour C->B")
    elif start[2] != end[2]:
        if end[2] == capacities[2]:
            moves.append("fill C")
        elif end[2] == 0:
            moves.append("empty C")
        else:
            if start[0] != end[0]:
                moves.append("pour A->C")
            else:
                moves.append("pour B->C")

print(moves)