import copy

def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return block

def dfs(stacks, goal, moves, visited):
    # Convert stacks to a tuple of tuples for immutability and hashing
    state_tuple = tuple(tuple(stack) for stack in stacks)
    
    # Check if the current state matches the goal state
    if stacks == goal:
        return moves
    
    # Avoid revisiting the same state
    if state_tuple in visited:
        return None
    visited.add(state_tuple)
    
    for source in range(3):
        if not stacks[source]:
            continue
        for destination in range(3):
            if source != destination:
                # Make a move
                block = move_block(stacks, source, destination)
                moves.append(f"Move {block} from {source + 1} to {destination + 1}")
                
                # Recurse
                result = dfs(copy.deepcopy(stacks), goal, moves, visited)
                if result:
                    return result
                
                # Undo the move
                move_block(stacks, destination, source)
                moves.pop()
    
    return None

initial_stacks = [['H', 'D'], ['E', 'J', 'A', 'C', 'F'], ['B', 'I', 'G', 'K']]
goal_stacks = [['B', 'E', 'I', 'J', 'K'], ['A', 'C', 'D', 'F'], ['G', 'H']]

moves = dfs(copy.deepcopy(initial_stacks), goal_stacks, [], set())
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")