import sys

# In case the user runs this script in a different environment,
# redirect stdout to a file to avoid mixed language output.
# For this platform, this is not strictly necessary as the output is captured.
# original_stdout = sys.stdout
# sys.stdout = open('output.txt', 'w', encoding='utf-8')

def solve():
    """
    Solves the problem by deriving and analyzing a recurrence relation for S(n).
    """
    p = 23627
    N = 510
    num_A = 203

    print("Step 1 & 2: Derive the recurrence relation for S(n).")
    print("Let S(n) be the number of valid 2 x n colorings.")
    print("Let N be the total number of colors, N = 510.")
    print("Let |A| be the number of restricted colors, |A| = 203.")
    print("A 2x3 rectangle is monochromatic with a restricted color c if columns i, i+1, i+2 are all (c,c).")
    print("By setting up a system of linear recurrences, we can find a recurrence for S(n):")
    print("S(n) = N^2 * S(n-1) - S(n-2) + (N^2 - |A|) * S(n-3)")
    print("-" * 20)

    print(f"Step 3: Simplify the recurrence modulo p = {p}.")
    N2_mod_p = (N * N) % p
    A_mod_p = num_A % p
    print(f"N^2 = {N*N}")
    print(f"N^2 mod {p} = {N2_mod_p}")
    print(f"|A| mod {p} = {A_mod_p}")

    if N2_mod_p == A_mod_p:
        print("We notice that N^2 = |A| (mod p). Let's call this value alpha.")
        alpha = N2_mod_p
        print(f"alpha = {alpha}")
        print("The recurrence for s_k = S(k) mod p becomes:")
        print("s_k = alpha * s_{k-1} - s_{k-2} + (alpha - alpha) * s_{k-3}")
        print("s_k = alpha * s_{k-1} - s_{k-2} for k >= 3.")
    else:
        print("N^2 and |A| are not equivalent mod p, the problem would be much harder.")
        # The logic below depends on this equivalence.
        return

    print("-" * 20)
    
    print("Step 4: Analyze the simplified recurrence.")
    print("We need the first few terms to analyze the sequence s_k = S(k) mod p.")
    s0 = 1 # S(0) is for an empty grid, which has one coloring (the empty one).
    s1 = N2_mod_p # S(1) = N^2
    s2 = (s1 * s1) % p # S(2) = (N^2)^2
    
    # Check if the recurrence s_k = alpha*s_{k-1} - s_{k-2} holds for k=2.
    # s2_check = (alpha * s1 - s0) % p
    # print(f"s2 from recurrence: {s2_check}, s2 actual: {s2}")
    # It does not hold for k=2, so the sequence s_k starts following the recurrence from k=3.
    # This means the sequence (s_k) for k>=1 has a state vector (s_k, s_{k+1}) that evolves
    # linearly for k>=2.
    
    print("The sequence s_k = S(k) mod p for k >= 1 follows the recurrence s_{k+2} = alpha * s_{k+1} - s_k for k>=1.")
    # My thought process showed s_k satisfies it for k>=3, or s_{k+2} for k>=1.
    s3_check = (alpha * s2 - s1) % p
    print(f"Let's verify for k=3: s3 = alpha * s2 - s1 = {alpha} * {s2} - {s1} = {s3_check}")
    S3_actual = (pow(N, 6, p * p) - num_A) # Use a large modulus to avoid issues
    s3_actual = S3_actual % p
    print(f"Actual s3 = (N^2)^3 - |A| mod p = ({pow(N2_mod_p,3,p)} - {A_mod_p}) mod p = {s3_actual}")
    print("The check passes. The sequence s_k for k>=1 satisfies s_{k+2} = alpha*s_{k+1}-s_k for k>=1.")
    print("This sequence is periodic. The period T_p divides p+1.")
    p_plus_1 = p + 1
    print(f"The period divides p+1 = {p_plus_1}.")
    print("-" * 20)

    print("Step 5: Calculate n modulo the period.")
    # n = 23626 * (23628^100 - 23628^50)
    # p+1 = 23628
    # n = (p-1) * ((p+1)^100 - (p+1)^50)
    # n mod (p+1) = (p-1 mod p+1) * (0 - 0) mod (p+1) = 0
    print(f"The given n is 23626 * (23628^100 - 23628^50).")
    print(f"Let x = p+1 = 23628. Then n = (x-2) * (x^100 - x^50).")
    print(f"n = x^101 - x^51 - 2*x^100 + 2*x^50.")
    print(f"Since all exponents are >= 1, n is a multiple of x = p+1 = {p_plus_1}.")
    print(f"So, n mod {p_plus_1} = 0.")
    print("Since the period T_p divides p+1, n is also a multiple of T_p.")
    print("-" * 20)
    
    print("Step 6: Compute the final answer.")
    print(f"We need to compute s_n. Since n is a multiple of the period T_p, s_n = s_{{T_p}} (for n large enough).")
    print("From the recurrence s_{k+2} = alpha * s_{k+1} - s_k for k>=1, we have s_2 = alpha * s_1 - s_{T_p} due to periodicity.")
    print(f"Therefore, s_{{T_p}} = alpha * s_1 - s_2.")
    final_eq_rhs_val = (alpha * s1 - s2) % p
    
    print("Let's plug in the numbers:")
    print(f"s_1 = {s1}")
    print(f"s_2 = {s2}")
    print(f"alpha = {alpha}")
    print("\nThe final equation is:")
    print(f"S(n) mod {p} = (alpha * s_1 - s_2) mod {p}")
    print(f"S(n) mod {p} = ({alpha} * {s1} - {s2}) mod {p}")
    term1 = (alpha * s1) % p
    print(f"S(n) mod {p} = ({term1} - {s2}) mod {p}")
    print(f"S(n) mod {p} = {final_eq_rhs_val}")
    
    print("\n" + "="*20)
    print("Final answer computation:")
    final_answer = final_eq_rhs_val
    print(f"The value S(23626*(23628^100-23628^50)) mod 23627 is {final_answer}.")
    print("="*20)

solve()

# sys.stdout = original_stdout
