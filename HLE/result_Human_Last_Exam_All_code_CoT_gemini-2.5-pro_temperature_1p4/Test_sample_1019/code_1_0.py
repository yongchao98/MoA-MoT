import math

def calculate_hamiltonicity_threshold():
    """
    Calculates the d-threshold for Hamiltonicity for H_n U G(n, p).

    The problem defines a graph H_n with minimum degree d >= n/2 - eta,
    for 1/2 <= eta <= n/64. The d-threshold is the minimum probability p
    such that H_n U G(n, p) is asymptotically almost surely Hamiltonian.

    The formula for the threshold is p = (2*eta - 1) / (n * (n - 2*eta)).

    To provide a concrete example, we choose n=1000.
    For this n, the valid range for eta is [0.5, 15.625]. We choose eta = 10.
    """
    n = 1000
    eta = 10

    # Verify that the chosen values are within the specified range
    if not (0.5 <= eta <= n / 64):
        print(f"Error: The chosen eta={eta} is not in the valid range for n={n}.")
        return

    # Numerator of the formula
    numerator = 2 * eta - 1
    
    # Denominator of the formula
    denominator_part1 = n
    denominator_part2 = n - 2 * eta
    denominator = denominator_part1 * denominator_part2
    
    # Calculate the threshold probability p
    p = numerator / denominator

    # Output the result in the format of the full equation
    print("The formula for the d-threshold for Hamiltonicity is:")
    print("p = (2 * eta - 1) / (n * (n - 2 * eta))")
    print("\nFor n = {} and eta = {}:".format(n, eta))
    
    # Print the equation with substituted values
    print("p = (2 * {} - 1) / ({} * ({} - 2 * {}))".format(eta, n, n, eta))
    print("p = {} / ({} * {})".format(numerator, denominator_part1, denominator_part2))
    print("p = {} / {}".format(numerator, denominator))
    print("p =", p)

if __name__ == '__main__':
    calculate_hamiltonicity_threshold()
