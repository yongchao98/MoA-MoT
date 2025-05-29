def get_moves_from_states(states, capacities):
    moves = []
    for i in range(1, len(states)):
        prev_state = states[i - 1]
        curr_state = states[i]
        
        # Determine the move
        if prev_state[0] != curr_state[0]:
            if curr_state[0] == capacities[0]:
                moves.append("fill A")
            elif curr_state[0] == 0:
                moves.append("empty A")
            else:
                if prev_state[1] != curr_state[1]:
                    moves.append("pour B->A")
                else:
                    moves.append("pour C->A")
        
        elif prev_state[1] != curr_state[1]:
            if curr_state[1] == capacities[1]:
                moves.append("fill B")
            elif curr_state[1] == 0:
                moves.append("empty B")
            else:
                if prev_state[0] != curr_state[0]:
                    moves.append("pour A->B")
                else:
                    moves.append("pour C->B")
        
        elif prev_state[2] != curr_state[2]:
            if curr_state[2] == capacities[2]:
                moves.append("fill C")
            elif curr_state[2] == 0:
                moves.append("empty C")
            else:
                if prev_state[0] != curr_state[0]:
                    moves.append("pour A->C")
                else:
                    moves.append("pour B->C")
    
    return moves

capacities = (4, 11, 11)
solution_states = [(0, 11, 0), (0, 11, 11), (4, 7, 11), (0, 7, 11), (4, 3, 11), (0, 3, 11), (3, 0, 11), (4, 0, 10), (0, 0, 10), (4, 0, 6), (0, 0, 6), (4, 0, 2)]
moves = get_moves_from_states(solution_states, capacities)
print(moves)