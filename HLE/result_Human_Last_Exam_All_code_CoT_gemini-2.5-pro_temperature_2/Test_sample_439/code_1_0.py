import itertools
import hashlib

def run_program_in_P(program_string):
    """
    Simulates running a program from the primitive recursive language P.
    Every program in P is guaranteed to halt and produce an output.
    To simulate this, we use a hash function to deterministically map
    the program string to an integer output. This ensures that for a given
    'program', the 'output' is always the same.
    """
    # The actual process of parsing and running a P-program would be here.
    # We use a hash as a stand-in for a deterministic, halting computation.
    # We take the integer value of the hash and get its remainder
    # when divided by 1000 to keep the outputs in a reasonable range.
    hash_object = hashlib.sha256(program_string.encode('utf-8'))
    hex_digest = hash_object.hexdigest()
    return int(hex_digest, 16) % 1000

def find_K(n):
    """
    Computes K(n) by brute-forcing all possible programs in increasing order of length.
    """
    # Let's assume our 'programs' are made of simple alphanumeric characters.
    # This represents the set of valid symbols in language P.
    charset = 'abcdefghijklmnopqrstuvwxyz0123456789'
    
    print(f"Searching for the shortest program that outputs n = {n}...\n")
    
    # Iterate through possible program lengths, starting from 1.
    for length in itertools.count(1):
        print(f"Checking all programs of length {length}...")
        # Generate all possible programs of the current length.
        # This is equivalent to `itertools.product(charset, repeat=length)`.
        for p_tuple in itertools.product(charset, repeat=length):
            program = "".join(p_tuple)
            
            # Since every program in P halts, this function call is guaranteed to return.
            output = run_program_in_P(program)
            
            if output == n:
                print("\n--- Found a solution! ---")
                print(f"The shortest program that outputs {n} is '{program}'.")
                print(f"The length of this program is {length}.")
                # Print the final equation as requested.
                print(f"K({n}) = {length}")
                return length

# --- Main Execution ---
if __name__ == "__main__":
    # The number we want to find the shortest program for.
    target_n = 42
    
    # Since the algorithm is computable, this function will halt and return a value.
    find_K(target_n)
