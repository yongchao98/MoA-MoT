import math

def solve():
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra system.

    The number of equilibria in the system is determined by the number of possible
    subsets of species that can coexist. As derived in the explanation, by varying
    the system parameters, any subset of the N species can form an equilibrium.
    This includes the trivial equilibrium where all species are extinct.

    The total number of subsets of a set with N elements is 2^N.
    Therefore, the system can have 2^N possible equilibria.
    """
    try:
        n_str = input("Enter the number of species (N): ")
        N = int(n_str)

        if N < 0:
            print("The number of species must be a non-negative integer.")
            return

        # Calculate the number of equilibria
        num_equilibria = 2**N

        print(f"For a system with N = {N} species, the number of possible equilibria is 2^N.")
        # The final equation is 2^N = result. We print each part.
        print(f"The calculation is: 2^{N} = {num_equilibria}")
        
        # The final answer format is not specified for a formula, so we will output the formula as a string.
        # If N is a specific number, the final answer would be the result of 2^N.
        # Since N is a variable, the answer is the formula itself.
        # Let's assume the user wants the formula as the answer.
        # For example, if N=3, the answer is 8. If N is not given, the answer is 2^N.
        # The prompt implies N is a variable, so the answer is a formula.
        # Let's output the formula as the final answer.
        # print("<<<2**N>>>") # This would be ambiguous.
        # Let's output the calculated number for the given N.
        # print(f"<<<{num_equilibria}>>>")

    except ValueError:
        print("Invalid input. Please enter an integer for N.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve()