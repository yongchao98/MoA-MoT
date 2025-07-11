def calculate_max_equilibria():
    """
    Calculates the maximum number of equilibria for the generalized Lotka-Volterra system.

    The maximum number of equilibria is 2^N, where N is the number of species.
    This script prompts the user for N and prints the result of the calculation.
    """
    try:
        # Prompt the user to enter the number of species, N.
        n_str = input("Please enter the number of species (N): ")
        n = int(n_str)

        if n < 0:
            print("Error: The number of species (N) cannot be negative.")
        else:
            # The maximum number of equilibria follows the equation 2^N.
            base = 2
            exponent = n
            result = base ** exponent

            print(f"\nFor a system with N = {n} species, the maximum number of possible equilibria is 2^N.")
            print("The calculation is shown below, with each number in the equation printed:")
            # Printing the equation in the format: base ** exponent = result
            print(f"{base} ** {exponent} = {result}")

    except ValueError:
        print("Invalid input. Please enter a valid integer for N.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    calculate_max_equilibria()