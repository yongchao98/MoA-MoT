import itertools

def find_shortest_program(target_n):
    """
    This function demonstrates the computability of K(n) for a primitive 
    recursive language P.

    Our toy language 'P' consists of programs made from the digits '0'-'9' and
    the '+' operator. For example, '3+5' is a program of length 3 that outputs 8.

    The function works by brute-force search:
    1. It tries all programs of length 1.
    2. Then all programs of length 2.
    3. And so on, until it finds the first program that outputs `target_n`.

    Because all programs in a primitive recursive language are guaranteed to
    halt, this search is guaranteed to be able to evaluate every program and
    eventually find a solution.
    """
    # Define the alphabet for our simple programming language P
    alphabet = '0123456789+'
    
    # K(n) is the length of the program. We start searching from length 1.
    program_length = 1
    
    print(f"Searching for the shortest program that outputs: {target_n}\n")
    
    while True:
        print(f"--- Checking all programs of length {program_length} ---")
        
        # Generate all possible strings of the current length from the alphabet
        possible_programs = itertools.product(alphabet, repeat=program_length)
        
        for p_tuple in possible_programs:
            program_str = "".join(p_tuple)
            
            # This is our 'interpreter' for language P. It must be safe and
            # always terminate. Using eval is safe here because our limited
            # alphabet prevents non-terminating code. We also catch syntax
            # errors for invalid programs like '5++' or '+'.
            try:
                # Security check: ensure the program only contains our alphabet chars
                # This is a safeguard if using a more powerful evaluator.
                if not all(c in alphabet for c in program_str):
                    continue

                result = eval(program_str)

                # Check if this program produces the target number
                if result == target_n:
                    print(f"\nSuccess! Found a shortest program.")
                    print(f"Program: '{program_str}'")
                    print(f"Output: {result}")
                    print(f"The value of K({target_n}) is the program length.")
                    # Final Answer formatting
                    print(f"{target_n} = K(n), length of \"{program_str}\" is {program_length}")
                    return

            except (SyntaxError, NameError, TypeError, ZeroDivisionError):
                # This program is not valid, so we ignore it.
                continue
                
        # If no program of the current length worked, try the next length
        program_length += 1

# Let's compute K(n) for a sample number.
# A small number like 21 is good, as it can be formed by '12+9' (len 4)
# or '20+1' (len 4) or '3*7' (if '*' were in the alphabet). With just '+',
# we expect something like '9+9+3' which is longer, or maybe '19+2'.
# Let's try to find K(21).
find_shortest_program(21)