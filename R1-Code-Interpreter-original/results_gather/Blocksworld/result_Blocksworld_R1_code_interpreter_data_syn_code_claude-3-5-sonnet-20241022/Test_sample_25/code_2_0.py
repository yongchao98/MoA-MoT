class BlocksWorld:
    def __init__(self):
        self.stacks = [
            ['A', 'G'],        # Stack 1
            ['F', 'E'],        # Stack 2
            ['C'],            # Stack 3
            ['H', 'J', 'I', 'B', 'D']  # Stack 4
        ]
        self.moves = []

    def print_state(self):
        for i, stack in enumerate(self.stacks):
            print(f"Stack {i+1}: {' '.join(stack)}")
        print()

    def get_top_block(self, stack_num):
        if not self.stacks[stack_num-1]:
            return None
        return self.stacks[stack_num-1][-1]

    def move(self, from_stack, to_stack):
        if not self.stacks[from_stack-1]:
            return False
        block = self.stacks[from_stack-1].pop()
        self.stacks[to_stack-1].append(block)
        move_str = f"Move {block} from {from_stack} to {to_stack}"
        self.moves.append(move_str)
        return True

    def solve(self):
        # First, get D to stack 3
        self.move(4, 3)  # Move D from 4 to 3
        
        # Get B to stack 1
        self.move(4, 1)  # Move B from 4 to 1
        
        # Get I to stack 3
        self.move(4, 3)  # Move I from 4 to 3
        
        # Get J to stack 2
        self.move(4, 2)  # Move J from 4 to 2
        
        # Get H to stack 3
        self.move(4, 3)  # Move H from 4 to 3
        
        # Get G to stack 2
        self.move(1, 2)  # Move G from 1 to 2
        
        # Get E is already in stack 2
        
        # Get F to stack 3
        self.move(2, 3)  # Move F from 2 to 3
        
        # Get C to stack 2
        self.move(3, 2)  # Move C from 2 to 3

        print('\n'.join(self.moves))

# Run the solution
world = BlocksWorld()
world.solve()