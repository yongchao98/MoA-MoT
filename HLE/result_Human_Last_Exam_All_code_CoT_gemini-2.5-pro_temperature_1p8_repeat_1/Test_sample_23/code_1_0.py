def solve():
    """
    Calculates the number of non-admissible integers for the given matrix problem.
    The user is prompted to enter two positive integers, a and b.
    """
    try:
        # Prompting the user for input
        a_str = input("Enter the value for positive integer a: ")
        b_str = input("Enter the value for positive integer b: ")
        a = int(a_str)
        b = int(b_str)

        if a <= 0 or b <= 0:
            print("The inputs a and b must be positive integers.")
            return

    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # To simplify the logic, we can assume a <= b without loss of generality.
    if a > b:
        a, b = b, a

    if a == 1:
        # Case 1: min(a,b) = 1
        num_inadmissible = a * b
        print(f"For a={a} and b={b}, since min(a, b) = 1, the number of non-admissible integers is a * b.")
        # Printing each number in the final equation
        print(f"Number = {a} * {b}")
        print(f"Result: {num_inadmissible}")
    else:
        # Case 2: a,b >= 2
        num_inadmissible = a + b - 2
        print(f"For a={a} and b={b}, since min(a, b) >= 2, the number of non-admissible integers is a + b - 2.")
        # Printing each number in the final equation
        print(f"Number = {a} + {b} - 2")
        print(f"Result: {num_inadmissible}")

# Execute the main function
solve()