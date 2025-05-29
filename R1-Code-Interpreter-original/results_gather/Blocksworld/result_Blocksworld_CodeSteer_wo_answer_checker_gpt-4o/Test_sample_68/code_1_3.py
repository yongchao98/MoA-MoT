def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return block

def dfs(stacks, goal, moves):
    if stacks == goal:
        return moves
    
    for source in range(3):
        if not stacks[source]:
            continue
        for destination in range(3):
            if source != destination:
                # Make a move
                block = move_block(stacks, source, destination)
                moves.append(f"Move {block} from {source + 1} to {destination + 1}")
                
                # Recurse
                result = dfs(stacks, goal, moves)
                if result:
                    return result
                
                # Undo the move
                move_block(stacks, destination, source)
                moves.pop()
    return None

initial_stacks = [['H', 'D'], ['E', 'J', 'A', 'C', 'F'], ['B', 'I', 'G', 'K']]
goal_stacks = [['B', 'E', 'I', 'J', 'K'], ['A', 'C', 'D', 'F'], ['G', 'H']]

moves = dfs(initial_stacks, goal_stacks, [])
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")