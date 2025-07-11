import math
from functools import reduce

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return abs(a)

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of integers."""
    if not numbers:
        return 0
    return reduce(gcd, numbers)

def generate_s_from_a(A):
    """Generates the set S from a given sequence A."""
    if not A:
        return set()
    n = len(A)
    generated_s = set()
    for i in range(n):
        current_gcd = A[i]
        generated_s.add(current_gcd)
        for j in range(i + 1, n):
            current_gcd = gcd(current_gcd, A[j])
            generated_s.add(current_gcd)
    return generated_s

def restore_a(S):
    """
    Attempts to restore a sequence A from a set S.
    S is provided as a list of integers.
    """
    if not S:
        print("Input set S is empty. No solution.")
        print("-1")
        return

    s_set = set(S)
    s_list = sorted(list(s_set))
    print(f"Given set S = {s_list}")

    # Step 1: Check the necessary condition.
    g = s_list[0]
    s_gcd = gcd_list(s_list)
    if g != s_gcd:
        print(f"Necessary condition failed: min(S)={g} is not equal to gcd(S)={s_gcd}.")
        print("-1")
        return
    
    print(f"Necessary condition passed: min(S) == gcd(S) == {g}.")

    # Step 2: Attempt a construction for A.
    # Handle the trivial case where S has one element.
    if len(s_list) == 1:
        A = s_list
        print(f"Found a trivial solution A = {A}")
        print("Final Answer Equation:")
        print(*A)
        return

    # Propose a simple construction: interleave g with other elements.
    s_prime = [s for s in s_list if s != g]
    candidate_A = []
    if s_prime:
        candidate_A.append(s_prime[0])
        for i in range(1, len(s_prime)):
            candidate_A.append(g)
            candidate_A.append(s_prime[i])
    
    if not candidate_A: # Should not happen if len > 1
        print("Construction failed.")
        print("-1")
        return
        
    # Step 3: Verify the constructed sequence.
    print(f"Trying a construction: A = {candidate_A}")
    generated_s = generate_s_from_a(candidate_A)

    if generated_s == s_set:
        print("Construction successful. A valid sequence is:")
        print("Final Answer Equation:")
        # The prompt requires printing each number of the equation
        output_str = " ".join(map(str, candidate_A))
        print(output_str)
    else:
        # A different construction might work. As shown in the analysis, sometimes none exists.
        # Another simple guess: A = S \ {g}
        candidate_A_2 = s_prime
        print(f"First construction failed. Generated S was {sorted(list(generated_s))}")
        print(f"Trying another construction: A = {candidate_A_2}")
        generated_s_2 = generate_s_from_a(candidate_A_2)
        if generated_s_2 == s_set:
            print("Construction successful. A valid sequence is:")
            print("Final Answer Equation:")
            output_str = " ".join(map(str, candidate_A_2))
            print(output_str)
        else:
            print(f"Second construction failed. Generated S was {sorted(list(generated_s_2))}")
            print("-1")

if __name__ == '__main__':
    # Example where a simple construction works.
    example_s = [6, 12, 18, 30]
    restore_a(example_s)
    
    print("\n" + "="*40 + "\n")
    
    # Example showing the necessary condition is not sufficient
    # (this algorithm will fail to find a solution).
    counter_example_s = [1, 2, 3, 5, 6, 10, 15]
    restore_a(counter_example_s)
