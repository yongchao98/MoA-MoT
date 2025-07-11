import sys

def solve():
    """
    Solves the Bit++ problem by reading statements from standard input,
    calculating the final value of x, and printing the step-by-step equation.
    """
    try:
        # Read the number of statements
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())

        x = 0
        # This list will store parts of the equation for the "show your work" output
        equation_parts = ["0"]

        for _ in range(n):
            # Read a single statement
            statement = sys.stdin.readline().strip()
            if not statement:
                continue
            
            # Check for the '+' character to determine the operation
            if '+' in statement:
                x += 1
                equation_parts.append("+ 1")
            else:
                x -= 1
                equation_parts.append("- 1")
        
        # This part fulfills the "output each number in the final equation" instruction
        # by printing the full calculation.
        equation_str = " ".join(equation_parts)
        print(f"Calculation: {equation_str} = {x}")
        
        # This part prints the final answer as required by the problem statement.
        print(f"Final Value: {x}")

    except (ValueError, IndexError):
        # Handle potential empty lines or malformed input gracefully
        print("Invalid input.", file=sys.stderr)

solve()