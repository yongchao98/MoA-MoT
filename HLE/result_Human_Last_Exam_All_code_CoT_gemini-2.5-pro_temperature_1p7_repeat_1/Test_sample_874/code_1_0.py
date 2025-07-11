import sys

# It's better to increase recursion limit for deep calculations if needed, though this solution is iterative.
# sys.setrecursionlimit(2000)

memo = {}

def f(a, b, c, d):
    """
    Calculates the length of the Ducci sequence for a given tuple.
    Includes the initial and final (0,0,0,0) tuples.
    Uses memoization to speed up calculations for repeated subproblems.
    """
    t = tuple(sorted((a,b,c,d)))
    if t in memo:
        return memo[t]

    original_tuple = (a, b, c, d)
    current_tuple = original_tuple
    # Using a set to detect cycles, although for Z_n this process always terminates.
    # The set of tuples is used to avoid recomputing f for permutations or cyclic shifts encountered during recursion.
    path = {tuple(sorted(current_tuple))}
    count = 1
    
    while True:
        if current_tuple == (0, 0, 0, 0):
            break
        
        current_tuple = (abs(current_tuple[0] - current_tuple[1]),
                         abs(current_tuple[1] - current_tuple[2]),
                         abs(current_tuple[2] - current_tuple[3]),
                         abs(current_tuple[3] - current_tuple[0]))
        count += 1
        
        # Check if subtuple is already computed
        sorted_sub_tuple = tuple(sorted(current_tuple))
        if sorted_sub_tuple in memo:
            count += memo[sorted_sub_tuple] - 1
            break

    # Memoize results for all intermediate sorted tuples on the path
    # to speed up future computations. This part is omitted for simplicity in final code
    # as the search space is narrowed down, but it's essential for a broad search.
    memo[t] = count
    return count


def solve():
    """
    Finds the tuple (a,b,c,d) and computes the required expression.
    """
    limit = 10_000_000

    # Step 1: Generate Tribonacci numbers
    trib = [0, 0, 1]
    while trib[-1] <= limit:
        trib.append(trib[-1] + trib[-2] + trib[-3])

    # Step 2: Find the maximum length M
    max_len = 0
    best_candidates = []

    # Iterate through primitive tuples of the form (t_i, t_{i+1}, t_{i+2}, 0)
    for i in range(len(trib) - 2):
        # We need t_{i+2} to be positive to form a non-trivial tuple
        if trib[i+2] == 0:
            continue
        # The largest element must be within the limit for scaling
        if trib[i+2] > limit:
            break

        a, b, c = trib[i], trib[i+1], trib[i+2]
        
        length = f(a, b, c, 0)
        
        if length > max_len:
            max_len = length
            best_candidates = [(i, length)]
        elif length == max_len:
            best_candidates.append((i, length))
    
    # Step 3: Find the candidate with the minimum sum
    min_sum = float('inf')
    final_tuple_primitive_index = -1

    for i, length in best_candidates:
        max_val = trib[i+2]
        if max_val == 0: continue

        s = limit // max_val
        current_sum = s * (trib[i] + trib[i+1] + trib[i+2])

        if current_sum < min_sum:
            min_sum = current_sum
            final_tuple_primitive_index = i

    # Step 4: Determine the final tuple and calculate the result
    i = final_tuple_primitive_index
    max_val = trib[i+2]
    s = limit // max_val

    # We choose the lexicographically smallest permutation of the scaled tuple.
    # The elements are {0, s*t_i, s*t_{i+1}, s*t_{i+2}}.
    # The smallest permutation is (0, s*t_i, s*t_{i+1}, s*t_{i+2}).
    a = 0
    b = s * trib[i]
    c = s * trib[i+1]
    d = s * trib[i+2]
    
    print(f"Maximum length M found: {max_len}")
    print(f"Best primitive tuple found is based on index i={i}, which is ({trib[i]}, {trib[i+1]}, {trib[i+2]}, 0)")
    print(f"Scaling factor used: {s}")
    print(f"The tuple (a, b, c, d) with max length and minimal sum is chosen as the lexicographically smallest permutation:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"d = {d}")
    print(f"This tuple has sum a+b+c+d = {a+b+c+d}")

    # Step 5: Compute the expression a + b - c - d mod 1000
    result = (a + b - c - d) % 1000
    print("\nFinal calculation:")
    print(f"{a} + {b} - {c} - {d} mod 1000 = {result}")
    
    # Final answer in required format
    print(f"\n<<<{result}>>>")

solve()
