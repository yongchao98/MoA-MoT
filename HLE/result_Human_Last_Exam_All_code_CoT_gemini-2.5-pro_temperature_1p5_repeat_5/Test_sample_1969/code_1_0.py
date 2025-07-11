import sys

def solve():
    """
    Calculates the fair division of money c_k for n players using Shapley values.
    The script takes n as a command-line argument, or defaults to 4.
    """
    try:
        # Default to n=4 if no command-line argument is given.
        n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
        if n <= 1:
            print("Error: n must be an integer greater than 1.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Please provide a valid integer for n.", file=sys.stderr)
        return

    print(f"Calculating fair division for n = {n} players.")
    print("-" * 40)

    # M is the sum of the first n integers
    # Using integer division // as the result is always an integer
    M = n * (n + 1) // 2
    
    # S2 is the sum of the first n squared integers
    S2 = n * (n + 1) * (2 * n + 1) // 6

    # Total earnings of the grand coalition
    total_earnings = M**4
    
    print(f"Total earnings to be divided: ({M})^4 = {total_earnings}")
    print("-" * 40)
    print("The derived formula for the amount c_k for player p_k is:")
    print(f"c_k = k * M * (M^2 + S2 - M * k)")
    print(f"where M = n(n+1)/2 = {M}")
    print(f"and S2 = n(n+1)(2n+1)/6 = {S2}")
    print("-" * 40)

    calculated_sum = 0
    for k in range(1, n + 1):
        # We use the formula c_k = k * M * (M^2 + S2 - M*k).
        # All calculations use integer arithmetic to maintain precision.
        m_squared = M**2
        inner_term = m_squared + S2 - M * k
        c_k = k * M * inner_term
        
        calculated_sum += c_k
        
        print(f"For player p_{k}:")
        # As requested, output each number in the final equation.
        print(f"c_{k} = {k} * {M} * ({m_squared} + {S2} - {M} * {k})")
        print(f"c_{k} = {k} * {M} * ({inner_term})")
        print(f"c_{k} = {c_k}\n")

    print("-" * 40)
    print(f"Sum of all individual payments: {calculated_sum}")
    if calculated_sum == total_earnings:
        print("The sum of payments correctly equals the total earnings.")
    else:
        # This case should not be reached with the correct formula.
        print(f"Error: Sum of payments ({calculated_sum}) does not match total earnings ({total_earnings}).")

if __name__ == '__main__':
    solve()