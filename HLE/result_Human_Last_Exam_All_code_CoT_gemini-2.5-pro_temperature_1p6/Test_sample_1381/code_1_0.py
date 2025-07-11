def solve():
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra dynamics.
    
    The user is prompted to enter the number of species (N).
    The code then calculates 2**N and prints the result.
    """
    try:
        n_str = input("Enter the number of species (N): ")
        n = int(n_str)
        if n < 0:
            print("Number of species must be a non-negative integer.")
            return

        num_equilibria = 2**n

        # The final equation is of the form: 2^N = result
        # We output each part of this equation.
        print(f"For N = {n} species, the number of possible equilibria is 2**{n} = {num_equilibria}.")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve()