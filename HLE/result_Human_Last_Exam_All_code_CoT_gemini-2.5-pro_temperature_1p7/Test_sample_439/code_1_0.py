import itertools

def find_shortest_program(target_n):
    """
    Simulates the computation of K(n) for a toy primitive recursive language.

    The language is defined as follows:
    - Programs are represented by nested tuples (Abstract Syntax Trees).
    - A constant value is a program: ('val', number). Length is 1.
      Our language's alphabet only includes primitive constants for 0, 1, and 2.
    - Addition of two programs is a new program: ('add', p1, p2). 
      Its length is 1 + length(p1) + length(p2).

    Since all programs are just combinations of constants and addition, they are
    guaranteed to halt, making K(n) computable for this language.
    """
    
    if not isinstance(target_n, int) or target_n < 0:
        print("Error: Input must be a non-negative integer.")
        return

    memo_eval = {}
    def evaluate(program):
        """Evaluates a program AST, with memoization for efficiency."""
        if program in memo_eval:
            return memo_eval[program]
        
        op = program[0]
        if op == 'val':
            result = program[1]
        elif op == 'add':
            result = evaluate(program[1]) + evaluate(program[2])
        else:
            # This case should not be reached with valid programs
            raise TypeError("Unknown program operation")
        
        memo_eval[program] = result
        return result

    def program_to_string(program):
        """Converts a program AST to a human-readable string equation."""
        op = program[0]
        if op == 'val':
            return str(program[1])
        elif op == 'add':
            return f"({program_to_string(program[1])} + {program_to_string(program[2])})"

    # A dictionary to store sets of programs, keyed by their length.
    programs_by_length = {}
    
    print(f"Searching for the shortest program that outputs {target_n}...")

    # Main loop: Iterate through all possible program lengths, starting from 1.
    for length in itertools.count(1):
        print(f"\n--- Generating programs of length {length} ---")
        
        current_length_programs = set()

        # Rule 1: Generate constant value programs.
        if length == 1:
            # Our toy language only has primitive constants 0, 1, and 2.
            for i in [0, 1, 2]:
                current_length_programs.add(('val', i))
        
        # Rule 2: Generate addition programs by combining shorter programs.
        # An 'add' program has length 1 + len(p1) + len(p2).
        # We iterate through all possible length combinations for p1 and p2.
        for len1 in range(1, length - 1):
            len2 = length - 1 - len1
            if len1 in programs_by_length and len2 in programs_by_length:
                for p1 in programs_by_length[len1]:
                    for p2 in programs_by_length[len2]:
                        # To avoid duplicates like (p1+p2) and (p2+p1) producing
                        # the same structure, we enforce an order.
                        if id(p1) <= id(p2):
                           current_length_programs.add(('add', p1, p2))
                        
        if not current_length_programs:
            print("No new programs can be constructed for this length.")
            continue
        
        programs_by_length[length] = current_length_programs
        print(f"Found {len(programs_by_length[length])} unique program(s).")
        
        # Check if any of the newly generated programs produce the target number.
        # We sort for deterministic output, making the example repeatable.
        for p in sorted(list(programs_by_length[length]), key=str):
            result = evaluate(p)
            if result == target_n:
                print("\n==============================================")
                print(">>> Success! Found a shortest program. <<<")
                print(f"Target Number (n): {target_n}")
                print(f"Program Length K({target_n}): {length}")
                # The final equation, showing each number involved
                print(f"Final Equation: {program_to_string(p)} = {result}")
                print("==============================================")
                return

if __name__ == '__main__':
    # We will find the shortest program for the number 9 as an example.
    target_number = 9
    find_shortest_program(target_number)
