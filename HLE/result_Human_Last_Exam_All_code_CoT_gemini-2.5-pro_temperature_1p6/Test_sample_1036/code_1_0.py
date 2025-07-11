def solve():
    """
    This function explains the reasoning and prints the number of positive integers n
    that satisfy the given property.

    Let n be a positive integer such that n <= lcm(1, 2, ..., 100).
    The property is that n gives different remainders when divided by each of k = 2, 3, ..., 100.
    Let r_k = n % k. The set R = {r_2, r_3, ..., r_100} must contain 99 distinct values.

    By definition, 0 <= r_k < k.
    For any k in {2, ..., 100}, we have r_k < k <= 100.
    So all 99 distinct remainders r_k must belong to the set {0, 1, ..., 99}.
    This means the set of remainders R must be {0, 1, ..., 99} with exactly one element missing.
    Let's call the missing element j. So, R = {0, 1, ..., 99} \ {j}.

    Case 1: j = 0.
    The set of remainders R is {1, 2, ..., 99}.
    For each k in {2, ..., 100}, the remainder r_k = n % k must be one of the values in R and must satisfy r_k < k.
    This forces a unique assignment:
    - r_2 must be < 2. The only choice from R is 1. So r_2 = 1.
    - r_3 must be < 3 and different from r_2. The only choice from R \ {1} is 2. So r_3 = 2.
    - ...
    - r_k = k - 1 for all k in {2, ..., 100}.
    This gives the system of congruences n % k = k - 1, which is n = -1 (mod k) for k=2..100.
    By CRT, n = -1 (mod lcm(2, ..., 100)).
    Let L = lcm(1, ..., 100). The only solution in the range 1 <= n <= L is n = L - 1.
    This gives one solution.

    Case 2: j > 0.
    A contradiction arises. For example, consider j = 1.
    The set of remainders R is {0, 2, 3, ..., 99}.
    Let's analyze the properties of n.
    - For k > j+1, r_k is forced to be k-1. So n = -1 (mod k) for k > j+1.
    - For j=1, n = -1 (mod k) for k in {3, ..., 100}.
    - This means n+1 is a multiple of lcm(3, ..., 100). Since 4 is in {3, ..., 100},
      lcm(3, ..., 100) is even. So n+1 is even, which means n is odd.
    - If n is odd, then r_2 = n % 2 = 1.
    - But for j=1, the set of remainders is {0, 2, ..., 99}. The value 1 is not in this set.
    - So, r_2 cannot be 1. This is a contradiction. There are no solutions for j=1.
    
    Similar contradictions can be derived for all other j > 0.
    
    Thus, there is only one integer n that satisfies the given conditions.
    The number of such integers is 1.
    """
    
    # The final equation is simply the count.
    number_of_integers = 1
    print(f"The number of such positive integers is: {number_of_integers}")

solve()