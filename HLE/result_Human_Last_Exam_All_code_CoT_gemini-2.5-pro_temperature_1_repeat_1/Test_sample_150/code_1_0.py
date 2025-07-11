import sys

def solve():
    """
    This function solves the Bit++ problem by reading a series of statements
    and calculating the final value of a variable 'x', initially 0.
    It also constructs and prints the full equation representing the operations.
    """
    try:
        # Read the first line to get the number of statements.
        # Handle the case of no input.
        line = sys.stdin.readline()
        if not line:
            print("0 = 0")
            return
        n = int(line)
    except (IOError, ValueError):
        # Handle empty input or non-integer input gracefully
        print("0 = 0")
        return

    x = 0
    # The equation starts with the initial value of x.
    equation_parts = ["0"]

    # Loop n times to process each statement
    for _ in range(n):
        # Read one statement from input and remove leading/trailing whitespace
        statement = sys.stdin.readline().strip()

        # If the statement is empty, skip it.
        if not statement:
            continue

        # Check for the presence of the '+' character to identify an increment.
        if '+' in statement:
            x += 1
            equation_parts.append("+")
            equation_parts.append("1")
        # Otherwise, it must be a decrement operation.
        else:
            x -= 1
            equation_parts.append("-")
            equation_parts.append("1")

    # Join all parts of the equation with spaces for readability.
    equation_str = " ".join(equation_parts)

    # Print the final equation and the result.
    print(f"{equation_str} = {x}")

# Execute the solution
solve()