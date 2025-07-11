import sys

def solve():
    """
    Calculates the number of possible equilibria for the generalized Lotka-Volterra system with N species.

    The number of possible equilibria corresponds to the number of subsets of species that can coexist,
    which has been shown to be 2^N.
    """
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py N")
        print("Please provide the number of species N as a command-line argument.")
        sys.exit(1)

    try:
        N = int(sys.argv[1])
        if N < 0:
            print("Error: The number of species (N) must be a non-negative integer.")
            sys.exit(1)
    except ValueError:
        print("Error: Please provide a valid integer for N.")
        sys.exit(1)

    # The number of possible equilibria is 2^N.
    num_equilibria = 2**N

    # As requested, printing the components of the final equation.
    print(f"For N = {N}, the number of possible equilibria is:")
    print(f"2^{N} = {num_equilibria}")

if __name__ == "__main__":
    solve()
