import math

# Memoization cache to store results for fractions we've already calculated
memo = {}

def count_hackenbush_pieces(m, d):
    """
    Recursively calculates the number of pieces for the canonical game
    representing the fraction m/d in Red-Blue Hackenbush.
    """
    # Use a tuple (m, d) as the key for the memoization cache
    original_m, original_d = m, d
    
    # Always simplify the fraction first
    common_divisor = math.gcd(m, d)
    m //= common_divisor
    d //= common_divisor

    if (m, d) in memo:
        return memo[(m, d)]

    # --- Base Cases ---
    # Case 1: The fraction is an integer (denominator is 1)
    if d == 1:
        # For an integer n, the game is n red pieces. Number of pieces is n.
        print(f"N({m}/{d}) = {abs(m)} (Integer)")
        memo[(m, d)] = abs(m)
        return abs(m)
    
    # Case 2: The fraction is 0
    if m == 0:
        # The game for 0 is the empty game, with 0 pieces.
        print(f"N(0/{d}) = 0 (Zero)")
        memo[(m, d)] = 0
        return 0

    # --- Recursive Step ---
    # Separate the integer and fractional parts
    integer_part = m // d
    fractional_m = m % d

    if fractional_m == 0:
        # This case should be handled by the d=1 base case after simplification
        # but is included for completeness.
        num_pieces = abs(integer_part)
        memo[(m, d)] = num_pieces
        return num_pieces

    # For a fraction y = m'/d', the rule is N(y) = 1 + N((m'-1)/d') + N((m'+1)/d')
    # First, calculate pieces for the integer part (if any)
    if integer_part != 0:
        print(f"Decomposing N({original_m}/{original_d}): N({integer_part}) + N({fractional_m}/{d})")
        pieces_for_integer = count_hackenbush_pieces(integer_part, 1)
    else:
        pieces_for_integer = 0
        
    # Now, calculate pieces for the fractional part
    left_m, left_d = fractional_m - 1, d
    right_m, right_d = fractional_m + 1, d

    # Get simplified versions for printing
    s_left_m, s_left_d = left_m // math.gcd(left_m, left_d), left_d // math.gcd(left_m, left_d)
    s_right_m, s_right_d = right_m // math.gcd(right_m, right_d), right_d // math.gcd(right_m, right_d)
    
    print(f"N({fractional_m}/{d}) = 1 + N({left_m}/{left_d}) + N({right_m}/{right_d})")
    print(f"         = 1 + N({s_left_m}/{s_left_d}) + N({s_right_m}/{s_right_d})")

    # Recursively call for the left and right options
    pieces_for_left = count_hackenbush_pieces(left_m, left_d)
    pieces_for_right = count_hackenbush_pieces(right_m, right_d)

    # The number of pieces for the fractional part is 1 (the root piece)
    # plus the pieces for the left and right options.
    pieces_for_fraction = 1 + pieces_for_left + pieces_for_right
    print(f"N({fractional_m}/{d}) = 1 + {pieces_for_left} + {pieces_for_right} = {pieces_for_fraction}")
    
    # Total pieces is sum of integer and fractional parts
    total_pieces = pieces_for_integer + pieces_for_fraction
    
    # Store result in cache before returning
    memo[(m, d)] = total_pieces
    
    return total_pieces

if __name__ == "__main__":
    numerator = 13
    denominator = 16
    print(f"Calculating the number of pieces for the fraction {numerator}/{denominator}:")
    print("-" * 30)
    final_answer = count_hackenbush_pieces(numerator, denominator)
    print("-" * 30)
    print(f"The total number of pieces needed to represent {numerator}/{denominator} is: {final_answer}")
    print(f"\nThe final equation for the top level is:")
    
    # To satisfy the "final equation" requirement, let's re-calculate the direct children
    n1 = 12
    d1 = 16
    n2 = 14
    d2 = 16
    
    # Ensure they are in the cache from the previous run
    val1 = memo[(n1 // math.gcd(n1,d1), d1 // math.gcd(n1,d1))]
    val2 = memo[(n2 // math.gcd(n2,d2), d2 // math.gcd(n2,d2))]
    
    print(f"N({numerator}/{denominator}) = 1 + N({n1}/{d1}) + N({n2}/{d2})")
    print(f"           = 1 + N(3/4) + N(7/8)")
    print(f"           = 1 + {val1} + {val2}")
    print(f"           = {final_answer}")