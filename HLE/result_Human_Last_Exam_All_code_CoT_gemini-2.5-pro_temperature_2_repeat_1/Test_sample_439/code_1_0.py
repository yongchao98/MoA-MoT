import itertools

def is_valid_p_program(program_str):
    """
    Checks if a string is a valid program in our toy language.
    A valid program must start with 'Z' and be followed only by 'S's.
    """
    if not program_str or program_str[0] != 'Z':
        return False
    for char in program_str[1:]:
        if char != 'S':
            return False
    return True

def execute_p_program(program_str):
    """
    Executes a valid program from our toy language.
    The output value is simply the number of 'S' characters.
    """
    # Assumes program_str is valid.
    return program_str.count('S')

def compute_k(n):
    """
    Computes K(n) for our toy language by brute-forcing through
    all possible programs, ordered by length.
    """
    # The alphabet of our language P
    alphabet = ['Z', 'S']
    
    # Iterate through all possible program lengths, starting from 1.
    L = 1
    while True:
        print(f"Searching for programs of length {L}...")
        
        # Generate all possible strings of length L from the alphabet.
        program_candidates = itertools.product(alphabet, repeat=L)
        
        for candidate_tuple in program_candidates:
            program_str = "".join(candidate_tuple)
            
            # Check if the generated string is a syntactically valid program.
            if is_valid_p_program(program_str):
                print(f"  - Testing valid program: '{program_str}'")
                
                # Execute the program. This is guaranteed to halt.
                output = execute_p_program(program_str)
                print(f"    - Program output: {output}")
                
                # If the output matches our target n, we've found the shortest program.
                if output == n:
                    print(f"\nSUCCESS: Found the shortest program for n = {n}.")
                    print(f"The program is '{program_str}' with length {L}.")
                    # Return the length L, which is the value of K(n).
                    return L
        
        # If no program of length L worked, try the next length.
        L += 1

# --- Main execution ---
# Let's compute K(n) for a sample number n.
# For example, to output 3, we need the program "ZSSS", which has length 4.
# So, we expect K(3) = 4.
target_n = 3
k_of_n = compute_k(target_n)

print("\n--- Final Result ---")
print(f"The equation is: K({target_n}) = {k_of_n}")
