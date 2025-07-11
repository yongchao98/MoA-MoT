import itertools

def compute_k_for_toy_language(target_n):
    """
    Computes K(n) for a simple toy language.

    K(n) is the length of the shortest program that outputs n.
    The language has two operations:
    - 'S': successor (adds 1)
    - 'D': double (multiplies by 2)
    A program is a sequence of these operations applied to a base value of 0.
    """
    if target_n == 0:
        print("A program of length 1, 'D', can produce 0.")
        print("Equation: D(0) = 0")
        print("K(0) = 1")
        return

    length = 1
    while True:
        # Generate all possible programs of the current length
        # 'SD' means apply S first, then D. D(S(0))
        programs = itertools.product('SD', repeat=length)

        # print(f"Checking programs of length {length}...")
        for program_tuple in programs:
            program_str = "".join(program_tuple)

            # Evaluate the program
            value = 0
            # Build a string to represent the equation for printing
            expr_str = "0"
            for op in program_str:
                if op == 'S':
                    value = value + 1
                    expr_str = f"S({expr_str})"
                elif op == 'D':
                    value = value * 2
                    expr_str = f"D({expr_str})"

            # Check if we found the target number
            if value == target_n:
                print(f"Found the shortest program for n={target_n}: '{program_str}'")
                print(f"The final equation is: {expr_str} = {value}")
                print(f"The length of this program is {length}.")
                return

        # If no program of this length worked, try the next length
        length += 1

# --- Main execution ---
# Set the integer 'n' for which we want to find K(n)
n = 10

print(f"Searching for the shortest program to output n = {n}...")
compute_k_for_toy_language(n)