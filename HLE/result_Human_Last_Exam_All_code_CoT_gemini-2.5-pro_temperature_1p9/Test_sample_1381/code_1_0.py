def calculate_equilibria():
    """
    Calculates the number of possible equilibria for the given generalized Lotka-Volterra system.
    """
    try:
        # Prompt the user to enter the number of species, N.
        n_str = input("Enter the number of species (N): ")
        n = int(n_str)
        if n < 0:
            print("Error: Number of species (N) cannot be negative.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a non-negative integer for N.")
        return

    # The number of equilibria is 2^N.
    base = 2
    exponent = n
    
    # Calculate the result.
    # We use a simple loop for exponentiation to handle large numbers gracefully
    # and to clearly demonstrate the process, though `base ** exponent` is more direct.
    num_equilibria = 1
    for _ in range(exponent):
        num_equilibria *= base

    # Print the final result and the numbers used in the equation.
    print(f"\nFor a system with N = {n} species, the maximum number of possible equilibria is {base}^{exponent} = {num_equilibria}.")
    print("The numbers in the final equation are:")
    print(f"Base: {base}")
    print(f"Exponent (N): {exponent}")
    print(f"Result: {num_equilibria}")


if __name__ == "__main__":
    calculate_equilibria()
