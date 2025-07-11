import collections

def find_shortest_program(n: int):
    """
    Simulates finding the shortest program for a number n in a toy language.

    This function demonstrates that K(n) is computable by performing a
    brute-force search for the shortest program.

    Our toy language starts with the number 1 and can apply two operations
    that are guaranteed to terminate:
    - Multiply by 2
    - Add 3

    A 'program' is a sequence of these operations. The 'length' of the
    program is the number of operations. We use a Breadth-First Search (BFS)
    to find the shortest sequence of operations to generate n, which is
    guaranteed to find the shortest path because it explores the search
    space layer by layer.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: Input must be an integer greater than or equal to 1.")
        return

    if n == 1:
        print("The shortest program for n=1 is the starting value itself.")
        print("Program length K(1) = 0")
        print("Final Equation: 1")
        return

    # A queue for the BFS. Each item is a tuple:
    # (current_number, string_representing_the_equation, length_of_program)
    queue = collections.deque([(1, "1", 0)])
    
    # A set to store visited numbers to prevent cycles and redundant work.
    visited = {1}

    while queue:
        current_num, equation_str, length = queue.popleft()

        # Define the possible operations (our "language primitives")
        operations = [
            {'op': '*', 'val': 2, 'str': '* 2'},
            {'op': '+', 'val': 3, 'str': '+ 3'}
        ]

        for op_info in operations:
            if op_info['op'] == '*':
                next_num = current_num * op_info['val']
            else: # op is '+'
                next_num = current_num + op_info['val']

            # Construct the string for the new equation
            next_equation_str = f"({equation_str} {op_info['str']})"
            next_length = length + 1
            
            # Check if we have found the target number
            if next_num == n:
                print(f"Found the shortest program for n={n}.")
                print(f"Program length K({n}) = {next_length}")
                # The final equation shows the 'program' and all numbers used.
                print(f"Final Equation: {next_equation_str} = {n}")
                return

            # If the number has not been visited, add it to the queue and visited set.
            # We add a pruning condition to avoid exploring huge numbers unnecessarily.
            if next_num not in visited and next_num < n * 2:
                visited.add(next_num)
                queue.append((next_num, next_equation_str, next_length))
    
    print(f"Could not find a program to generate n={n} with the given operations.")

# Example: Find the shortest program for a target number.
# Let's find the shortest program for n = 37.
target_number = 37
find_shortest_program(target_number)
