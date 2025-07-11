import itertools

def compute_K_n(target_n):
    """
    Demonstrates the computability of K(n) for a primitive recursive language.

    Our toy language 'P' consists of strings that are valid Python arithmetic
    expressions using integers, '+', and '*'. Every program in this language
    halts (evaluation either succeeds or raises a syntax error).

    K(n) is the length of the shortest program string that evaluates to n.
    """
    print(f"Searching for the shortest program that outputs: {target_n}\n")
    
    # The alphabet of our toy programming language P
    alphabet = "0123456789+*"
    
    # Iterate through possible program lengths, starting from 1
    for length in itertools.count(start=1):
        print(f"--- Checking programs of length {length} ---")
        
        # Generate all possible programs of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for p_tuple in possible_programs:
            program_str = "".join(p_tuple)
            
            # A valid program cannot start or end with an operator,
            # and cannot contain leading zeros in numbers (e.g., '05').
            if program_str.startswith(('+', '*')) or program_str.endswith(('+', '*')):
                continue
            if any(p in program_str for p in ['+0', '*0']) and len(program_str) > 1:
                if program_str not in ['0+0', '0*0']: # edge cases
                    # This is a simplification to avoid multi-digit numbers with leading zeros
                    # A full parser would be needed for perfect handling.
                    if program_str[program_str.find('0')-1] in '+*':
                         continue
            
            try:
                # Run the program (evaluate the expression)
                # This is safe because all programs are guaranteed to halt.
                result = eval(program_str)
                
                # Check if the output matches our target
                if result == target_n:
                    print(f"\nSUCCESS: Found the shortest program for {target_n}.")
                    print(f"Program: '{program_str}'")
                    print(f"Length (K({target_n})): {length}")
                    
                    # To satisfy the "output each number in the final equation" requirement
                    # We re-format the program string with spaces for readability
                    equation = program_str.replace("+", " + ").replace("*", " * ")
                    print(f"Final Equation: {equation} = {result}")
                    return
                    
            except (SyntaxError, NameError, TypeError):
                # This string is not a valid program in our language, so we ignore it.
                continue
                
# Let's find K(n) for a sample number.
# A small number is chosen so the search completes quickly.
target_number = 42
compute_K_n(target_number)
