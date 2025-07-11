import itertools
import re

def compute_k_simulation(target_n, charset, max_len=10):
    """
    Simulates the computation of K(n) for a given n in a toy language.
    The language consists of expressions formed by the given character set.
    """
    print(f"Searching for the shortest program to output {target_n}...")
    print(f"Using the language alphabet: {list(charset)}\n")

    # 1. Iterate through possible program lengths, starting from 1.
    for length in range(1, max_len + 1):
        print(f"Checking all {len(charset)**length} programs of length {length}...")
        
        # 2. Generate all possible programs of the current length.
        program_generator = itertools.product(charset, repeat=length)

        for p in program_generator:
            program_str = "".join(p)

            # Heuristic to skip obviously invalid programs before calling eval()
            if program_str.startswith(('*', '+')) or program_str.endswith(('*', '+')):
                continue
            if '**' in program_str or '++' in program_str or '*+' in program_str or '+*' in program_str:
                continue

            try:
                # 3. Run the program. This is guaranteed to halt in our language.
                result = eval(program_str)

                # 4. Check if the output matches the target.
                if result == target_n:
                    print("\n--- Found ---")
                    
                    # Split the program string to format the equation nicely
                    parts = re.split(r'([*+])', program_str)
                    equation = " ".join(parts) + f" = {result}"
                    
                    print(f"Shortest program found: '{program_str}'")
                    print(f"Final Equation: {equation}")
                    print(f"The length of this shortest program is {length}.")
                    print(f"Therefore, for this toy language, K({target_n}) = {length}.")
                    return

            except (SyntaxError, NameError, TypeError):
                # This program was not valid syntax, so we ignore it.
                continue
    
    print(f"\nCould not find a program to output {target_n} within length {max_len}.")

# --- Main execution ---
# Define our toy language's alphabet (its primitive components)
LANGUAGE_CHARSET = '4*+' 
# Define the target number we want to produce
TARGET_NUMBER = 64

# Run the simulation
compute_k_simulation(TARGET_NUMBER, LANGUAGE_CHARSET)