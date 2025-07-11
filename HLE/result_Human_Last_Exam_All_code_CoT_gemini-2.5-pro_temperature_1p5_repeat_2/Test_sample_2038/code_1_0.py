import math

def solve():
    """
    Calculates the number of 2-bridge knots with crossing number <= 13
    that admit two disjoint non-parallel minimal genus Seifert surfaces.
    """
    
    print("The problem is to count 2-bridge knots K(p,q) up to crossing number 13 with a special property.")
    print("This property is equivalent to the knot's continued fraction being of the form [2a_1, ..., 2a_n] with a_i >= 1 and n >= 2.")
    print("For the link to be a knot, n must be even.")
    print("The crossing number C = 2 * (a_1 + ... + a_n). So C must be even, and C >= 4.")
    print("Let N = C/2 = a_1 + ... + a_n. We need to count compositions of N into an even number of parts, up to reversal.")
    print("A knot and its mirror are counted as one.")
    print("-" * 20)

    total_knots = 0
    knot_counts_per_c = []

    # Crossing number C must be an even integer. Smallest is 2*1+2*1=4.
    # We check C from 4 to 12.
    for C in range(4, 14, 2):
        N = C // 2
        
        # A_e_N is the number of compositions of N into an even number of parts.
        # Formula: 2**(N-2) for N >= 2.
        if N < 2:
            A_e_N = 0
        else:
            A_e_N = 2**(N - 2)

        # P_e_N is the number of palindromic compositions of N into an even number of parts.
        # Formula: 2**(N/2 - 1) if N is even and N >= 2, else 0.
        P_e_N = 0
        if N >= 2 and N % 2 == 0:
            P_e_N = 2**(N // 2 - 1)
            
        # Total distinct knots for crossing number C.
        # This counts pairs of (composition, its_reverse) as one object.
        # This corresponds to counting (chiral knot, its_mirror) as one,
        # and amphichiral knots as one.
        num_knots_C = (A_e_N + P_e_N) // 2
        
        print(f"For crossing number C = {C}:")
        print(f"  N = C/2 = {N}")
        print(f"  Number of compositions of {N} into an even number of parts (A_e) = {A_e_N}")
        print(f"  Number of palindromic such compositions (P_e) = {P_e_N}")
        print(f"  Number of knots = (A_e + P_e) / 2 = ({A_e_N} + {P_e_N}) / 2 = {num_knots_C}")
        print("-" * 20)
        
        total_knots += num_knots_C
        knot_counts_per_c.append(str(num_knots_C))

    # Format the final output string
    equation = " + ".join(knot_counts_per_c)
    print("Total number of such knots is the sum over the allowed crossing numbers:")
    print(f"{equation} = {total_knots}")

solve()
print("<<<19>>>")