import itertools

def run_simple_P_program(program_code):
    """
    Runs a program from our simple primitive recursive language 'P'.
    - A valid program starts with 'Z' and is followed only by 'S's.
    - 'Z' outputs 0. Each 'S' acts as a successor function, adding 1.
    - Raises a ValueError for invalid program syntax.
    """
    if not program_code:
        raise ValueError("Program cannot be empty.")
    if program_code[0] != 'Z':
        raise ValueError("Program must start with 'Z'.")
    
    value = 0
    for i, char in enumerate(program_code[1:]):
        if char != 'S':
            raise ValueError(f"Invalid character '{char}' at position {i+1}.")
        value += 1
    return value

def compute_K(n):
    """
    Computes K(n) by searching for the shortest program that outputs n.
    This is a computable process because run_simple_P_program always halts.
    """
    if not isinstance(n, int) or n < 0:
        print("Input must be a non-negative integer.")
        return

    print(f"Computing K({n})...")
    
    program_length = 1
    while True:
        print(f"\n--- Searching for programs of length {program_length} ---")
        # Generate all possible programs of the current length from the alphabet {'Z', 'S'}
        alphabet = ['Z', 'S']
        possible_programs = itertools.product(alphabet, repeat=program_length)
        
        for program_tuple in possible_programs:
            program = "".join(program_tuple)
            try:
                output = run_simple_P_program(program)
                print(f"Trying program '{program}'... It halts and outputs: {output}")
                if output == n:
                    print("\n--- Found a shortest program! ---")
                    # The final "equation" showing the result
                    print(f"K({n}) = {program_length}")
                    return program_length
            except ValueError as e:
                # This program has invalid syntax in our language 'P'
                print(f"Trying program '{program}'... It is not a valid program. ({e})")
        
        program_length += 1

# Let's compute K(n) for a sample value, e.g., n=2.
# The code will demonstrate the search process that proves its computability.
compute_K(2)
