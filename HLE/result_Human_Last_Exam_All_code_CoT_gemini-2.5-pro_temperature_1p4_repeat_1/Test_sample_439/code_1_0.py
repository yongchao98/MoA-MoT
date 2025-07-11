import sys

def find_shortest_program(target_n):
    """
    This function simulates the computation of K(n) for a toy language.
    K(n) is the length of the shortest program that outputs n.

    Our toy language is defined in the `programs` dictionary.
    'S(x)' can be thought of as a successor function (adds 1).
    'a*b' can be thought of as multiplication.
    The length of a program is its string length.
    """
    # A dictionary representing our primitive recursive language `P`.
    # It maps program strings to their integer outputs.
    # In a real scenario, we would generate and execute these programs.
    # Here, we pre-compute them for this demonstration.
    programs = {
        "0": 0,
        "S(0)": 1,
        "S(S(0))": 2,
        "S(S(S(0)))": 3,
        "2*2": 4,
        "S(S(S(S(0))))": 4,
        "2+3": 5,
        "S(S(S(S(S(0)))))": 5,
        "2*3": 6,
        "S(S(S(S(S(S(0))))))": 6,
        "3+4": 7,
        "S(S(S(S(S(S(S(0)))))))": 7,
        "4*2": 8,
        "S(S(S(S(S(S(S(S(0))))))))": 8,
        "3*3": 9,
        "S(S(S(S(S(S(S(S(S(0)))))))))": 9,
    }

    # We determine the maximum possible length to search.
    # In a true implementation, this loop would be unbounded,
    # but we know it will halt.
    max_len = 0
    for p in programs:
        if len(p) > max_len:
            max_len = len(p)

    # The brute-force algorithm:
    # Iterate through every possible program length, from 1 upwards.
    for length in range(1, max_len + 2):
        shortest_program_found = ""
        
        # Check all programs of the current length.
        for program_str, output_val in programs.items():
            if len(program_str) == length:
                # If a program produces the target number, we've found the shortest one
                # for this length. Since we are iterating length upwards, this is the
                # overall shortest.
                if output_val == target_n:
                    shortest_program_found = program_str
                    break # Stop checking other programs of this length
        
        if shortest_program_found:
            k_n = len(shortest_program_found)
            # This "equation" shows the result: K(n) = length
            print(f"To find K({target_n}):")
            print(f"- Searched programs of lengths 1 through {length}.")
            print(f"- Found shortest program: '{shortest_program_found}'")
            print(f"- The output of this program is {target_n}.")
            print(f"- The length of this program is {k_n}.")
            print("-" * 20)
            print(f"Result: K({target_n}) = {k_n}")
            return # Halt the entire function

    print(f"Could not find a program that outputs {target_n} in our defined language.")


# --- Main Execution ---
# Set the target number 'n' for which we want to compute K(n).
# Let's use n=9, because our language has a short program "3*3" (length 3)
# and a long one "S(S(S(S(S(S(S(S(S(0)))))))))" (length 25).
# The algorithm should correctly find the shorter one.
target_number = 9
find_shortest_program(target_number)
