class BlocksWorld:
    def __init__(self):
        # Initialize current state
        self.current_state = {
            1: ['D'],
            2: ['A', 'J'],
            3: ['H', 'F', 'C', 'E', 'I'],
            4: ['B', 'G']
        }
        
        # Goal state
        self.goal_state = {
            1: ['F', 'I'],
            2: ['A', 'B', 'C', 'D', 'G'],
            3: ['E', 'H', 'J']
        }
        
        self.moves = []
    
    def get_top_block(self, stack_num):
        if self.current_state[stack_num]:
            return self.current_state[stack_num][-1]
        return None
    
    def is_valid_move(self, from_stack, to_stack):
        return len(self.current_state[from_stack]) > 0
    
    def make_move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.current_state[from_stack].pop()
            self.current_state[to_stack].append(block)
            self.moves.append(f"Move {block} from {from_stack} to {to_stack}")
    
    def solve(self):
        # Strategy: First clear space by moving blocks to stack 4 if needed
        # Then build target stacks one by one
        
        # First, let's build stack 1 (F, I)
        while self.get_top_block(3) != 'I':
            self.make_move(3, 4)
        self.make_move(3, 1)  # Move I to stack 1
        
        while self.get_top_block(3) != 'F':
            self.make_move(3, 4)
        self.make_move(3, 1)  # Move F to stack 1
        
        # Now build stack 2 (A, B, C, D, G)
        # First get A from stack 2
        while self.get_top_block(2) != 'A':
            self.make_move(2, 3)
        
        # Now build the rest
        blocks_to_stack2 = ['B', 'C', 'D', 'G']
        for block in blocks_to_stack2:
            # Find block and move obstacles
            for stack in [1, 2, 3, 4]:
                while self.current_state[stack] and self.get_top_block(stack) != block and len(self.current_state[stack]) > 0:
                    self.make_move(stack, 3)
                if self.current_state[stack] and self.get_top_block(stack) == block:
                    self.make_move(stack, 2)
                    break
        
        # Finally, build stack 3 (E, H, J)
        blocks_to_stack3 = ['E', 'H', 'J']
        for block in blocks_to_stack3:
            for stack in [1, 2, 3, 4]:
                while self.current_state[stack] and self.get_top_block(stack) != block and len(self.current_state[stack]) > 0:
                    if stack != 3:
                        self.make_move(stack, 4)
                if self.current_state[stack] and self.get_top_block(stack) == block:
                    self.make_move(stack, 3)
                    break

# Create and solve
solver = BlocksWorld()
solver.solve()
print("<<<")
print("\n".join(solver.moves))
print(">>>")