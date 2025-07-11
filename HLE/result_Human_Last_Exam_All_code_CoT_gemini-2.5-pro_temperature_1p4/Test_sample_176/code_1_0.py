import sys

def solve_xpp_interpreter():
    """
    Reads an X++ program from standard input, executes it,
    and prints the calculation as an equation.
    """
    try:
        # Read the number of statements from the first line.
        num_statements_str = sys.stdin.readline()
        if not num_statements_str:
            return
        num_statements = int(num_statements_str)

        x_value = 0
        # Start building the equation string with the initial value of X.
        equation = ["0"]

        for _ in range(num_statements):
            statement = sys.stdin.readline().strip()
            
            # Check for increment or decrement operation and update the value and equation.
            if "++" in statement:
                x_value += 1
                equation.append("+")
                equation.append("1")
            elif "--" in statement:
                x_value -= 1
                equation.append("-")
                equation.append("1")
        
        # Add the final result to the equation.
        equation.append("=")
        equation.append(str(x_value))
        
        # Print the full equation.
        print(" ".join(equation))

    except (ValueError, IndexError):
        # Handle cases with improper input format.
        print("Invalid input: Please provide a valid X++ program format.")

solve_xpp_interpreter()