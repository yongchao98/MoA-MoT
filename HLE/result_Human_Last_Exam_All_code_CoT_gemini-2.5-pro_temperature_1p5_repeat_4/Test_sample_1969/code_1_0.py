import math

def calculate_shapley_division(n):
    """
    Calculates the Shapley value for each player in the described game.

    Args:
        n (int): The total number of players. Must be greater than 1.

    Prints:
        The fair share c_k for each player k from 1 to n,
        along with the details of the calculation.
    """
    if not isinstance(n, int) or n <= 1:
        print("Please provide an integer n > 1.")
        return

    # Calculate S_1 and S_2, the sum of first powers and second powers of integers up to n
    s1 = n * (n + 1) // 2
    s2 = n * (n + 1) * (2 * n + 1) // 6

    # Pre-calculate the coefficients of the formula for c_k
    s1_squared = s1 * s1
    s1_cubed = s1_squared * s1
    
    coeff_k1 = s1_cubed + s1 * s2
    coeff_k2 = s1_squared

    print(f"For n = {n} players:")
    print(f"The total earnings are ({n}({n+1})/2)^4 = {s1**4}")
    print(f"The formula for the share c_k of player p_k is: c_k = ({s1}^3 + {s1}*{s2})*k - {s1}^2*k^2")
    print(f"c_k = {coeff_k1}*k - {coeff_k2}*k^2")
    print("-" * 20)

    total_payout = 0
    for k in range(1, n + 1):
        k_squared = k * k
        
        # Calculate the two terms of the formula
        term1_val = coeff_k1 * k
        term2_val = coeff_k2 * k_squared
        
        # Calculate the final share for player k
        c_k = term1_val - term2_val
        total_payout += c_k
        
        # Print the detailed equation for c_k
        print(f"c_{k} = {coeff_k1} * {k} - {coeff_k2} * {k_squared}")
        print(f"   = {term1_val} - {term2_val}")
        print(f"   = {c_k}\n")
        
    print("-" * 20)
    print(f"The sum of all shares is: {total_payout}")
    print(f"This matches the total earnings: {s1**4}")


if __name__ == '__main__':
    # You can change the value of n here to see the results for a different number of players.
    try:
        n_players = int(input("Enter the number of players (n > 1): "))
        calculate_shapley_division(n_players)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
