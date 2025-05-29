class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]

    def is_valid_move(self, from_stack, to_stack):
        # Check if source stack has blocks
        if not self.stacks[from_stack]:
            return False
        return True

    def make_move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.stacks[from_stack].pop()
            self.stacks[to_stack].append(block)
            return True
        return False

    def is_goal(self):
        goal = ['A', 'B', 'C', 'D', 'E']
        return len(self.stacks[2]) == 5 and self.stacks[2] == goal

def solve_blocksworld():
    # Initial state
    initial_state = State([['E', 'B'], ['A', 'D', 'C'], []])
    moves = []
    
    # Helper function to get current top blocks
    def get_top_blocks(state):
        return [stack[-1] if stack else None for stack in state.stacks]
    
    # Solution steps
    current_state = initial_state
    
    # Step 1: Clear C (top of stack2)
    current_state.make_move(1, 0)  # Move C from stack2 to stack1
    moves.append("Move C from 2 to 1")
    
    # Step 2: Clear D
    current_state.make_move(1, 2)  # Move D from stack2 to stack3
    moves.append("Move D from 2 to 3")
    
    # Step 3: Move A to stack3
    current_state.make_move(1, 2)  # Move A from stack2 to stack3
    moves.append("Move A from 2 to 3")
    
    # Step 4: Move C to stack2
    current_state.make_move(0, 1)  # Move C from stack1 to stack2
    moves.append("Move C from 1 to 2")
    
    # Step 5: Move B to stack3
    current_state.make_move(0, 2)  # Move B from stack1 to stack3
    moves.append("Move B from 1 to 3")
    
    # Step 6: Move C to stack3
    current_state.make_move(1, 2)  # Move C from stack2 to stack3
    moves.append("Move C from 2 to 3")
    
    # Step 7: Move E to stack2
    current_state.make_move(0, 1)  # Move E from stack1 to stack2
    moves.append("Move E from 1 to 2")
    
    # Step 8: Move D to stack3
    current_state.make_move(2, 1)  # Move D from stack3 to stack2
    moves.append("Move D from 3 to 2")
    
    # Step 9: Move E to stack3
    current_state.make_move(1, 2)  # Move E from stack2 to stack3
    moves.append("Move E from 2 to 3")
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()