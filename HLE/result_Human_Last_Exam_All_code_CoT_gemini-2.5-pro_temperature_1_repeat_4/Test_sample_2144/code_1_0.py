import math

def solve_for_x0():
    """
    This function solves for the position x0 by first finding the parameter p
    and then substituting it into the parametric equation for x.
    """
    # The problem reduces to solving the quartic equation for p: p^4 - 18p - 27 = 0.
    # We look for an integer root by testing the divisors of 27.
    p = None
    # Divisors of -27 are ±1, ±3, ±9, ±27
    for i in [1, -1, 3, -3, 9, -9, 27, -27]:
        if i**4 - 18*i - 27 == 0:
            p = i
            break

    if p is not None:
        print(f"The integer root of the equation p^4 - 18p - 27 = 0 is p = {p}.")
        
        # Now, we calculate x0 using the parametric equation:
        # x(p) = -3(p + 1) / sqrt(2p + 3)
        
        # Breaking down the calculation as requested
        val_p_plus_1 = p + 1
        val_2p_plus_3 = 2 * p + 3
        
        if val_2p_plus_3 < 0:
            print("Error: The value under the square root is negative.")
            return

        val_sqrt = math.sqrt(val_2p_plus_3)
        numerator = -3 * val_p_plus_1
        x0 = numerator / val_sqrt
        
        print("\nCalculating x0 using the value of p:")
        print(f"x0 = -3 * (p + 1) / sqrt(2 * p + 3)")
        print(f"x0 = -3 * ({p} + 1) / sqrt(2 * {p} + 3)")
        print(f"x0 = -3 * ({val_p_plus_1}) / sqrt({val_2p_plus_3})")
        print(f"x0 = {numerator} / {val_sqrt}")
        print(f"The position x0 is: {x0}")

    else:
        print("No integer root was found for the parameter p.")

solve_for_x0()