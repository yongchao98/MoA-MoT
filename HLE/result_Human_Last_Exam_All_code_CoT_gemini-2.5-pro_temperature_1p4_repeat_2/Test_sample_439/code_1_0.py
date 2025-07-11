import collections

def run_toy_program(p: str) -> int:
    """
    Interpreter for our simple toy language P.
    A program is a string of 'i', 'a', 'd'. It must start with 'i'.
    'i': initialize to 0
    'a': add 1
    'd': double
    This function is guaranteed to halt for any string input.
    """
    if not p or p[0] != 'i':
        # Invalid programs output -1 for simplicity
        return -1
    
    accumulator = 0
    # Process commands after 'i'
    for command in p[1:]:
        if command == 'a':
            accumulator += 1
        elif command == 'd':
            accumulator *= 2
        else:
            # Invalid command in program
            return -1
    return accumulator

def compute_K(n: int):
    """
    Computes K(n) for our toy language P by performing a breadth-first search.
    This demonstrates the computability of K(n) for a primitive recursive language.
    """
    if n < 0:
        print("This toy language can only produce non-negative integers.")
        return

    # A queue for the breadth-first search. We start with the shortest possible program.
    # The BFS naturally finds the shortest program first because we explore programs
    # layer by layer, where each layer corresponds to a program length.
    queue = collections.deque(["i"])
    
    # A set to keep track of programs we've already queued to avoid redundant computations.
    visited = {"i"}

    print(f"Searching for the shortest program to output n = {n}...")
    
    while queue:
        # Get the next program from the queue
        program = queue.popleft()
        
        # Run the program to get its output. This step is guaranteed to halt.
        output = run_toy_program(program)
        
        # Check if the program produces our target number n
        if output == n:
            length = len(program)
            print("\n--- Found a Solution! ---")
            print(f"The shortest program that outputs {n} is: '{program}'")
            # The final equation K(n) = length
            print(f"K({n}) = {length}")
            return length

        # If not, generate the next possible programs by adding 'a' or 'd'
        # and add them to the queue if we haven't seen them before.
        for op in ['a', 'd']:
            next_program = program + op
            if next_program not in visited:
                visited.add(next_program)
                queue.append(next_program)

if __name__ == '__main__':
    # We will compute K(n) for a sample number.
    # A small number is chosen so the search completes quickly.
    target_number = 20
    compute_K(target_number)
