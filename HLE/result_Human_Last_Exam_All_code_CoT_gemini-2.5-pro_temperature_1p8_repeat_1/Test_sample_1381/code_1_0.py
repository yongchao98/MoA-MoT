import math

def calculate_possible_equilibria(N):
    """
    This function calculates the total number of possible equilibria for a generalized
    Lotka-Volterra system with N species by demonstrating that any subset of species
    can form an equilibrium for some choice of parameters.

    The total number of equilibria is the sum of the number of possible subsets of
    surviving species of every size k, from 0 to N. This is given by the sum of
    binomial coefficients C(N, k).

    Args:
        N (int): The number of species in the system.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: The number of species N must be a non-negative integer.")
        return

    print(f"For a system with N = {N} species:")
    print("An equilibrium state is defined by its unique subset of surviving species.")
    print("We can find parameters such that any of the 2^N possible subsets of species can form a stable equilibrium.")
    print("\nThe number of ways to form an equilibrium with k surviving species is given by the binomial coefficient 'N choose k'.")
    print("To find the total number of possible equilibria, we sum these counts for all possible k from 0 to N:")

    # This list will store each term C(N, k) of the sum
    num_equilibria_per_size = []
    for k in range(N + 1):
        # Calculate C(N, k), which is the number of subsets of size k
        count = math.comb(N, k)
        num_equilibria_per_size.append(count)

    # Calculate the total, which is guaranteed to be 2**N
    total_equilibria = sum(num_equilibria_per_size)
    
    # Format the final equation string showing each number
    equation_str = " + ".join(map(str, num_equilibria_per_size))

    print("\nFinal Equation:")
    print(f"{equation_str} = {total_equilibria}")
    print(f"\nThus, there are {total_equilibria} possible equilibria as parameters vary.")


if __name__ == '__main__':
    # Set the number of species N for the calculation.
    # You can change this value to see the result for a different number of species.
    N = 10
    calculate_possible_equilibria(N)
    print("\nThe final answer for the total number of equilibria is 2^N.")
    final_answer = 2**N
    # The output format is specific, so we print the final calculation.
    # Here is an example of what would be returned for N=10:
    print(f"\nFor N={N}, the answer is {final_answer}")
    print(f'<<<2^{N}>>>')