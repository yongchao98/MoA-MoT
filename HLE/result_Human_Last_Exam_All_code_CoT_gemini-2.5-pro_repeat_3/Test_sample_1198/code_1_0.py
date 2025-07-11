import decimal
import math

# Set precision for Decimal calculations. 50 should be sufficient for this depth.
decimal.getcontext().prec = 50

def find_a(n, interval, modulus, max_depth):
    """
    Recursively search for a valid interval for 'a'.
    
    Args:
        n (int): The current power to check.
        interval (tuple): A tuple (L, R) representing the current interval [L, R) for 'a'.
        modulus (int): The modulus (2 or 3).
        max_depth (int): The maximum depth for the search.
        
    Returns:
        bool: True if a valid interval is found down to max_depth, False otherwise.
    """
    if n > max_depth:
        # We found a sequence of choices that works up to max_depth.
        # This means a candidate interval exists.
        return True

    L, R = interval
    
    # If the interval is vanishingly small, prune this search path.
    if R - L < decimal.Decimal('1e-40'):
        return False

    # Determine the range for k = floor(a^n)
    try:
        min_val_pow_n = L ** n
        max_val_pow_n = R ** n
    except decimal.InvalidOperation:
        # This can happen if numbers get too large, effectively a failure.
        return False

    k_min = min_val_pow_n.to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    # We need to find integers k in [min_val_pow_n, max_val_pow_n)
    # The list of candidate k goes from floor(L^n) to floor(R^n)
    k_max = max_val_pow_n.to_integral_value(rounding=decimal.ROUND_FLOOR)
    if k_max == max_val_pow_n:
        # R^n is an integer, so floor(a^n) must be strictly less than R^n
        k_max -= 1
    
    # Find all possible values of k that satisfy the congruence
    possible_k = []
    # Start check from the first integer >= min_val_pow_n
    current_k = min_val_pow_n.to_integral_value(rounding=decimal.ROUND_CEILING)
    
    while current_k <= k_max:
        if current_k % modulus == n % modulus:
            possible_k.append(current_k)
        current_k += 1

    if not possible_k:
        # No valid integer k found in the range, this path fails.
        return False

    # Try each valid k and recurse
    for k in possible_k:
        # New interval for 'a' is [k^(1/n), (k+1)^(1/n))
        # We need to intersect it with the current interval [L, R)
        k_dec = decimal.Decimal(k)
        # Using n-th root for decimals
        new_L = max(L, k_dec ** (decimal.Decimal(1)/n))
        new_R = min(R, (k_dec + 1) ** (decimal.Decimal(1)/n))

        if new_L < new_R:
            if find_a(n + 1, (new_L, new_R), modulus, max_depth):
                # A valid sub-path was found, so this path is potentially valid.
                return True
    
    # All choices of k from this state led to a dead end.
    return False

def solve():
    """
    Runs the search for both modulo 2 and modulo 3 and prints the final answer.
    """
    # Start with a wide interval for a > 1. Let's start with a=1.
    # For n=1, floor(a) must be 1 mod 2/3. Let's assume floor(a)=1.
    # So the initial interval is [1, 2).
    initial_interval = (decimal.Decimal(1), decimal.Decimal(2))
    
    # Check for modulo 2
    # A search depth of 20 is sufficient to show that the search space collapses.
    result_mod2 = find_a(1, initial_interval, 2, 20)
    
    # Check for modulo 3
    result_mod3 = find_a(1, initial_interval, 3, 20)
    
    answer1 = "Yes" if result_mod2 else "No"
    answer2 = "Yes" if result_mod3 else "No"
    
    print(f"{answer1},{answer2}")

if __name__ == '__main__':
    solve()