from collections import deque
import copy

class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def __str__(self):
        return str(self.stacks)
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    def __hash__(self):
        return hash(str(self))

def get_valid_moves(state):
    moves = []
    for i in range(len(state.stacks)):
        if state.stacks[i]:  # if stack is not empty
            block = state.stacks[i][-1]  # get top block
            for j in range(len(state.stacks)):
                if i != j:  # can't move to same stack
                    moves.append((block, i, j))
    return moves

def apply_move(state, move):
    block, source, dest = move
    new_stacks = copy.deepcopy(state.stacks)
    new_stacks[source].pop()
    new_stacks[dest].append(block)
    return State(new_stacks)

def bfs_search():
    initial = State([['B', 'G', 'F'], ['A', 'E', 'C'], ['D']])
    goal = State([['D', 'F', 'G'], ['E'], ['A', 'B', 'C']])
    
    queue = deque([(initial, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        if current_state == goal:
            return path
            
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        for move in get_valid_moves(current_state):
            next_state = apply_move(current_state, move)
            move_str = f"Move {move[0]} from {move[1]+1} to {move[2]+1}"
            queue.append((next_state, path + [move_str]))
    
    return None

solution = bfs_search()
if solution:
    print('\n'.join(solution))