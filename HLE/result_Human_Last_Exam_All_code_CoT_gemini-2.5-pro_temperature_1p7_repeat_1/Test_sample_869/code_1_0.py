import math

def solve_probability():
    """
    Calculates the probability P = F/S based on the problem description.
    """
    # N = 25 items total
    # T = 5 types of items
    # 5 copies of each type
    # 5 individuals, each gets 5 items

    # Step 1: Calculate the total number of possible distributions (S)
    # This is the number of distinct ways to arrange 25 items with 5 groups of 5 identical items.
    # S = 25! / (5! * 5! * 5! * 5! * 5!)
    try:
        s_numerator = math.factorial(25)
        s_denominator = math.factorial(5)**5
        S = s_numerator / s_denominator
    except OverflowError:
        # For very large numbers, we can use floating point representation directly.
        # This will be handled by the division below.
        pass

    # Step 2: Calculate the number of favorable distributions (F)
    # A favorable distribution requires each individual to have a unique "dominant" item type.
    # This restrictive condition is only met if each individual 'i' possesses all 5 items of a unique type 'pi(i)'.
    # The number of ways to assign these unique types to the 5 individuals is the number of permutations of 5.
    # F = 5!
    F = math.factorial(5)

    # Step 3: Calculate the probability P = F/S
    # We use the direct formula to maintain precision.
    # P = F / S = 5! / (25! / (5!^5)) = (5!^6) / 25!
    P = F / S

    # Output the result in the specified format
    print(f"The total number of ways to distribute the items is S.")
    print(f"The number of favorable distributions is F.")
    print(f"The probability is P = F / S.")
    print(f"")
    print(f"Calculation:")
    print(f"F = 5! = {F}")
    # S can be very large, print it in scientific notation if needed
    print(f"S = 25! / (5!)^5 = {S:e}")
    print(f"P = {F} / {S:.0f} = {P:e}")
    print(f"")
    # To show the final equation with all numbers
    s_val = int(S)
    f_val = int(F)
    print(f"Final Equation:")
    print(f"P = {f_val} / {s_val}")
    
    # We can also express P as (5!^6) / 25!
    final_p_num = math.factorial(5)**6
    final_p_den = math.factorial(25)
    print(f"P = (5!)^6 / 25! = {final_p_num} / {final_p_den}")


solve_probability()
<<<P = 120 / 623360743125120>>>