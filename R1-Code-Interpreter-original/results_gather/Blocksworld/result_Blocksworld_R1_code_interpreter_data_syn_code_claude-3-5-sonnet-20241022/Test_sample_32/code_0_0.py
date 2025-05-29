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

    def __str__(self):
        return str(self.stacks)

def find_solution():
    # Initial state
    initial = State([['G', 'D', 'B', 'C'], ['F', 'A', 'E'], []])
    # Goal state
    goal = [['A', 'B', 'D', 'E', 'G'], ['C'], ['F']]
    
    moves = []
    current = initial
    
    # Move C from stack1 to stack2
    current.make_move(0, 2)
    moves.append("Move C from 1 to 3")
    
    # Move E from stack2 to stack3
    current.make_move(1, 2)
    moves.append("Move E from 2 to 3")
    
    # Move B from stack1 to stack2
    current.make_move(0, 1)
    moves.append("Move B from 1 to 2")
    
    # Move D from stack1 to stack3
    current.make_move(0, 2)
    moves.append("Move D from 1 to 3")
    
    # Move G from stack1 to stack2
    current.make_move(0, 1)
    moves.append("Move G from 1 to 2")
    
    # Move D from stack3 to stack1
    current.make_move(2, 0)
    moves.append("Move D from 3 to 1")
    
    # Move E from stack3 to stack1
    current.make_move(2, 0)
    moves.append("Move E from 3 to 1")
    
    # Move C from stack3 to stack2
    current.make_move(2, 1)
    moves.append("Move C from 3 to 2")
    
    # Move G from stack2 to stack1
    current.make_move(1, 0)
    moves.append("Move G from 2 to 1")
    
    # Move B from stack2 to stack1
    current.make_move(1, 0)
    moves.append("Move B from 2 to 1")
    
    # Move A from stack2 to stack1
    current.make_move(1, 0)
    moves.append("Move A from 2 to 1")
    
    # Move F from stack2 to stack3
    current.make_move(1, 2)
    moves.append("Move F from 2 to 3")
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()