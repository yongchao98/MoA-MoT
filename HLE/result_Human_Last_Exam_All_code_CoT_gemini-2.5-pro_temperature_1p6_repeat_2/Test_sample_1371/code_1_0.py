import math

def calculate_arrangements():
    """
    This function calculates the number of ways to arrange the people at the circular table
    based on the given constraints.
    """

    # Pre-calculate factorials for convenience
    f10 = math.factorial(10)
    f5 = math.factorial(5)
    f4 = math.factorial(4)

    # Case 1: Classicists are NOT adjacent to the Scientist-Mathematician (SM) block.
    # The two ethicists are used as buffers.
    # Block arrangement: SM - E1 - C - E2 (or E2, E1 swapped -> 2 ways)
    # W_sm = 360 * 10!
    # W_c = 5!
    ok_arrangements_factor = 2 * 360 * f5
    ok_arrangements = ok_arrangements_factor * f10

    # Case 2: Classicists ARE adjacent to the SM block (via Cassie's rule).
    # This forms a new mega-block {SMC}, which is arranged with the two ethicists.
    # Number of ways to arrange {SMC}, E1, E2 in a circle is (3-1)! = 2.

    # Subcase 2a: Classicist block joins the Scientist-end of the SM block.
    # This requires a female scientist at the S-end and Cassie at the C-end.
    # W_sm_S_female = 216 * 10!
    # W_c_cassie = 4!
    cassie_S_factor = 2 * 216 * f4
    cassie_S_arrangements = cassie_S_factor * f10

    # Subcase 2b: Classicist block joins the Mathematician-end of the SM block.
    # This requires a female mathematician at the M-end and Cassie at the C-end.
    # W_sm_M_female = 240 * 10!
    # W_c_cassie = 4!
    cassie_M_factor = 2 * 240 * f4
    cassie_M_arrangements = cassie_M_factor * f10

    # Total arrangements is the sum of all cases
    total_factor = ok_arrangements_factor + cassie_S_factor + cassie_M_factor
    total_arrangements = total_factor * f10
    
    # Print the breakdown of the calculation
    print("The total number of arrangements is calculated by summing three cases:")
    print("1. Standard arrangement (Classicists not next to Scientists/Mathematicians):")
    print(f"   (2 block arrangements) * (360 * 10!) * 5! = {ok_arrangements_factor} * 10!")
    print("\n2. Cassie sits next to a female Scientist:")
    print(f"   (2 block arrangements) * (216 * 10!) * 4! = {cassie_S_factor} * 10!")
    print("\n3. Cassie sits next to a female Mathematician:")
    print(f"   (2 block arrangements) * (240 * 10!) * 4! = {cassie_M_factor} * 10!")
    
    print("\nSumming the factors:")
    print(f"({ok_arrangements_factor} + {cassie_S_factor} + {cassie_M_factor}) * 10! = {total_factor} * 10!")
    
    print("\nFinal equation and answer:")
    print(f"{total_factor} * {f10} = {total_arrangements}")


calculate_arrangements()
<<<108288 * 10!>>>