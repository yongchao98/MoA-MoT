import itertools
import re

def compute_k_for_n(target_n):
    """
    Computes K(n) for a given n in our simple primitive recursive language P.

    Language P:
    - Programs are strings.
    - The alphabet is '0123456789+*'.
    - Rule: Multi-digit numbers are not allowed (e.g., '10' is invalid).
    - The program's output is its evaluation as a Python expression.
    """
    print(f"Attempting to compute K({target_n}) for our sample language P...")
    print("Language P consists of expressions using single digits ('0'-'9') and operators ('+', '*').")
    print(f"K({target_n}) is the length of the shortest program (string) that evaluates to {target_n}.\n")

    alphabet = '0123456789+*'
    length = 1

    # Loop indefinitely, increasing program length until a solution is found.
    while True:
        print(f"Searching for programs of length {length}...")
        
        # Generate all possible program strings of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for p_tuple in possible_programs:
            program_str = "".join(p_tuple)

            # --- Filter for our Language P's rules ---
            # Rule: No multi-digit numbers (e.g. '11', '23')
            if re.search(r'\d\d', program_str):
                continue
            
            try:
                # --- Run the program and check its output ---
                result = eval(program_str)
                
                if result == target_n:
                    # Found the shortest program since we search by increasing length
                    print("\n--- Found Shortest Program ---")
                    print(f"Program String: '{program_str}'")
                    
                    # Output all the numbers in the final equation as requested
                    equation = program_str.replace('*', ' * ').replace('+', ' + ')
                    print(f"Final Equation: {equation} = {target_n}")
                    
                    k_n = len(program_str)
                    print(f"The length of this program is {k_n}.")
                    print(f"Therefore, for our language P, K({target_n}) = {k_n}")
                    return

            except (SyntaxError, ZeroDivisionError):
                # This string is not a valid expression in our language (e.g., '++' or '*9').
                # We simply ignore it and continue the search.
                continue
        
        length += 1
        # Safety break for this demonstration to prevent excessively long runs.
        if length > 10:
            print("\nSearch is taking too long for this demonstration. Halting.")
            break

if __name__ == '__main__':
    # We will compute K(42). In our language, this should be '6*7' or '7*6', which has a length of 3.
    # This is shorter than any other simple combination.
    target_number = 42
    compute_k_for_n(target_number)