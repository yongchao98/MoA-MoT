def get_moves():
    # Initial state and goal state
    initial = [['B', 'A'], ['C', 'F', 'G'], ['D', 'I', 'E', 'J', 'H']]
    goal = [['A', 'E', 'I'], ['B', 'C', 'J'], ['D', 'F', 'G', 'H']]
    
    moves = []
    current = initial
    
    # Helper function to get current state as string
    def state_str(state):
        return ','.join([''.join(stack) for stack in state])
    
    # Helper function to make a move
    def make_move(state, from_stack, to_stack):
        if not state[from_stack]:
            return None
        new_state = [stack[:] for stack in state]
        block = new_state[from_stack].pop()
        new_state[to_stack].append(block)
        return new_state, f"Move {block} from {from_stack + 1} to {to_stack + 1}"

    # First, clear the blocks that are not in their final position
    # Move H to correct position
    moves.append("Move H from 3 to 2")
    current[2].pop()  # Remove H
    current[1].append('H')
    
    # Move J
    moves.append("Move J from 3 to 1")
    current[2].pop()  # Remove J
    current[0].append('J')
    
    # Move E
    moves.append("Move E from 3 to 2")
    current[2].pop()  # Remove E
    current[1].append('E')
    
    # Move I
    moves.append("Move I from 3 to 1")
    current[2].pop()  # Remove I
    current[0].append('I')
    
    # Move G and F to stack 3
    moves.append("Move G from 2 to 3")
    current[1].pop()  # Remove G
    current[2].append('G')
    
    moves.append("Move F from 2 to 3")
    current[1].pop()  # Remove F
    current[2].append('F')
    
    # Move remaining blocks to their final positions
    moves.append("Move C from 2 to 1")
    moves.append("Move B from 1 to 2")
    moves.append("Move A from 1 to 3")
    moves.append("Move I from 1 to 3")
    moves.append("Move C from 1 to 2")
    moves.append("Move J from 1 to 2")
    moves.append("Move A from 3 to 1")
    moves.append("Move I from 3 to 1")
    moves.append("Move E from 2 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

get_moves()