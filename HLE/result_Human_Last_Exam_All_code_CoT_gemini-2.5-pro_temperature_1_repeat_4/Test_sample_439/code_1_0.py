import itertools

def run_primitive_program(program_str):
    """
    Runs a program in our simple 'S'/'Z' primitive recursive language.
    'Z' is the constant 0.
    'S' is the successor function (+1).
    The program is evaluated from right to left.
    Example: 'SSZ' -> S(S(Z)) -> S(S(0)) -> S(1) -> 2
    """
    if not program_str:
        return None # Invalid program

    # Find the starting value, which must be 'Z' in our simple language
    if program_str.endswith('Z'):
        value = 0
        # Apply successor 'S' functions
        for char in reversed(program_str[:-1]):
            if char == 'S':
                value += 1
            else:
                # Invalid program format for this simple interpreter
                return None
        return value
    else:
        # A valid program in our simple language must end with a base case 'Z'
        return None

def compute_K(target_n):
    """
    Finds the length of the shortest program that outputs target_n.
    It does this by performing a brute-force search, which is possible
    because all programs are guaranteed to halt.
    """
    print(f"Searching for the shortest program to output n = {target_n}...")

    # 1. Iterate through possible program lengths, starting from 1
    length = 1
    while True:
        print(f"\nChecking all programs of length {length}...")
        
        # 2. Generate all possible programs of the current length
        # Our language alphabet is {'S', 'Z'}
        alphabet = ['S', 'Z']
        possible_programs = itertools.product(alphabet, repeat=length)

        for p_tuple in possible_programs:
            program = "".join(p_tuple)
            
            # 3. Run the program and check its output
            result = run_primitive_program(program)
            # print(f"  Running '{program}' -> Output: {result}") # Uncomment for verbose output

            if result == target_n:
                # 4. If a match is found, we are done
                print(f"\nSuccess! Found shortest program.")
                print(f"Program '{program}' outputs {result}.")
                print(f"The length of the shortest program is {length}.")
                # This demonstrates the computability of K(n)
                # In our simple language, K(n) = n + 1
                # Final Equation: K({target_n}) = {length}
                print(f"Final Equation: K({target_n}) = {length}")
                return length
        
        # 5. If no program of this length worked, try the next length
        length += 1

if __name__ == '__main__':
    # Let's find K(n) for n = 4
    target_number = 4
    compute_K(target_number)