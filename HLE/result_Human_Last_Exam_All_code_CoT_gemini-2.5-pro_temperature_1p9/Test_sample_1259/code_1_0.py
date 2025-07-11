import math

def solve():
    """
    This function reads a set of integers S, determines if a generating
    sequence A can exist, and if so, prints one such sequence.
    """
    try:
        # Read the input set S. Assuming space-separated integers on one line.
        s_input = list(map(int, input().split()))
        s_set = sorted(list(set(s_input)))
    except (IOError, ValueError):
        print("Invalid input. Please provide a space-separated list of integers.")
        return

    if not s_set:
        print("Input set S cannot be empty.")
        return
        
    n = len(s_set)
    g = s_set[0] # The minimum element in the sorted set

    # Check the necessary condition: min(S) must be the gcd of all elements in S.
    # This is equivalent to checking if the minimum element divides all other elements.
    is_valid_s = True
    for x in s_set:
        if x % g != 0:
            is_valid_s = False
            break

    if not is_valid_s:
        # The necessary condition is not met, so no solution exists.
        print("-1")
        return

    # The condition is also sufficient. We can construct a solution.
    # One valid construction is to interleave the minimum element 'g'
    # between all other elements.
    result_a = []
    if n == 1:
        result_a.append(g)
    else:
        # Let's create the sequence [s'_1, g, s'_2, g, ..., g, s'_k]
        # where s' are the other elements.
        for x in s_set:
            if x == g:
                continue
            result_a.append(x)
            result_a.append(g)
        
        # The last appended 'g' might be extra if s_set[-1] is not g.
        # Let's fix it by removing the last element if it is g
        # A simpler way is to just print g, then for every other number print it and g.
        result_a_final = [g]
        for x in s_set:
            if x != g:
                result_a_final.append(x)
                result_a_final.append(g)
        
        # My analysis shows A = [s'_1, g, s'_2, ... g, s'_{k-1}] is a correct construction.
        # This yields a list of length 2 * |S'| - 1 = 2 * (n-1) - 1 = 2n-3
        final_A = []
        others = [x for x in s_set if x != g]
        final_A.append(g)
        for val in others:
          final_A.append(val)
        
        # A simple valid solution of length |S| can be just g followed by the other elements.
        # A = [g, s'_1, s'_2, ...]
        # Let's verify: gcd(g,s'_i)=g. gcd(s'_i, s'_j) = multiple of g. We can't guarantee this gcd is in S.
        
        # Let's use the proven interleaved method:
        # A = [g, s'_1, g, s'_2, ...] No this is wrong.
        # A = [s'_1, g, s'_2, g, ...] is the one that works.
        
        A_final = []
        if n > 1:
            others = [x for x in s_set if x != g]
            A_final.append(others[0])
            for i in range(1, len(others)):
                A_final.append(g)
                A_final.append(others[i])
        else:
            A_final.append(g)

    # Let's try another simple one. [s1, g, s2, g ...]. This also works.
    reconstructed_A = []
    for x in s_set:
        reconstructed_A.append(x)
        reconstructed_A.append(g)
    # This simple construction has duplicates of g. Let's make it more elegant.
    if g in reconstructed_A:
        reconstructed_A.pop() # Remove the final extra g

    # The prompt demands printing the equation. A = [n1, n2, ...].
    # Just print the elements space separated.
    # The simplest construction that works is just printing g and then every other element.
    # [s_1, g, s_2] -> gcd(s1,s2) is not checked.
    # so A = [g, s1, s2 ...] fails
    
    # A simple construction is `[g, s_1, g, s_2, ...]` This might fail if `gcd(s_i)` is not itself. No, elements are from S.
    # Okay let's stick to the interleaved construction which I proved.
    if n > 1:
        answer_A = []
        other_elements = s_set[1:]
        answer_A.append(other_elements[0])
        for i in range(1, len(other_elements)):
            answer_A.append(g)
            answer_A.append(other_elements[i])
        print(*answer_A)

    else: # n == 1
        print(s_set[0])

# To make this runnable, here is how you would use it.
# You can uncomment the below lines and run the script.
# Provide input like "6 12 18" in the console.
# if __name__ == '__main__':
#    solve()
# Example 1: Input "6 12 18" -> Output: "12 6 18"
# Example 2: Input "2 3 6" -> Output: "-1"
