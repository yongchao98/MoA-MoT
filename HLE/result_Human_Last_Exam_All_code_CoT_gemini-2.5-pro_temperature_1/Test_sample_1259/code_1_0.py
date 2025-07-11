import math

def solve():
    """
    Solves the problem of restoring sequence A from set S.
    The input is hardcoded for demonstration purposes, but in a real
    scenario, it would be read from standard input.
    """
    
    # Example input. Replace this with user input if needed.
    # For instance:
    # m = int(input())
    # s_input = list(map(int, input().split()))
    
    # Case 1: A solvable case
    s_input = [6, 12, 18]
    # Case 2: An unsolvable case
    # s_input = [6, 10, 15] 
    # Case 3: A simple case
    # s_input = {6, 12}

    # The problem defines S as a set. We should work with unique elements.
    s_unique = sorted(list(set(s_input)))
    m = len(s_unique)

    if m == 0:
        print(-1)
        return

    def gcd_list(numbers):
        if not numbers:
            return 0
        result = numbers[0]
        for i in range(1, len(numbers)):
            result = math.gcd(result, numbers[i])
        return result

    # The necessary and sufficient condition is that the minimum element
    # of the set must be equal to the GCD of all elements in the set.
    g_min = s_unique[0]
    g_all = gcd_list(s_unique)

    if g_min != g_all:
        print(-1)
    else:
        # If the condition holds, a solution can be constructed.
        # One valid construction is to interleave the minimum element (g_min)
        # with all other elements of the set.
        if m == 1:
            A = [g_min]
        else:
            A = []
            # s_unique is sorted, so g_min is at index 0.
            # The other elements are from index 1 onwards.
            other_elements = s_unique[1:]
            for element in other_elements:
                A.append(element)
                A.append(g_min)
        
        # The prompt asks to output the restored sequence.
        # The phrasing "output each number in the final equation" is interpreted
        # as printing the sequence A, which is the "solution to the equation".
        # We print the length of A, followed by its elements.
        print(len(A))
        print(*A)

# Execute the solution
solve()