import itertools

def evaluate_program(program_str):
    """
    Evaluates a program in our toy language {'S', '0'}.
    A valid program is a string of 'S's followed by a single '0'.
    For example, 'SS0' evaluates to 2.
    Returns the integer output or None if the program syntax is invalid.
    """
    s_count = program_str.count('S')
    zero_count = program_str.count('0')

    # A valid program must contain exactly one '0' and it must be at the end.
    if zero_count == 1 and program_str.endswith('0') and (s_count + zero_count == len(program_str)):
        return s_count
    else:
        return None

def compute_K(n):
    """
    This function simulates the computation of K(n) for our toy language.
    It finds the length of the shortest program that outputs n.
    """
    if not isinstance(n, int) or n < 0:
        print("This demonstration only works for non-negative integers.")
        return

    print(f"Attempting to compute K({n})...")
    
    alphabet = ['S', '0']
    length = 1

    while True:
        # Generate all possible programs of the current length
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for program_tuple in possible_programs:
            program_str = "".join(program_tuple)
            
            # Run the program and get its output
            result = evaluate_program(program_str)
            
            # Check if this program produces our target number n
            if result == n:
                print("\n--- Found Solution ---")
                print(f"The shortest program that outputs {n} is '{program_str}'.")
                print(f"The length of this program is {length}.")
                print(f"Therefore, the equation is: K({n}) = {length}")
                
                # As requested, printing each number in the final equation K(n)=L
                print("\nThe numbers in the final equation are:")
                print(n)
                print(length)
                return

        # If no program of the current length worked, try the next length
        length += 1

# --- Main execution ---
# Let's find the value of K(n) for a sample number.
target_number = 4
compute_K(target_number)
