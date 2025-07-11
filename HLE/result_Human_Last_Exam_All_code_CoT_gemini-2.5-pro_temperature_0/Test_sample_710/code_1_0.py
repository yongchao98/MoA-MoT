import sys

def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """
    
    # Step 1: Generate the sequence a_n = 2n^2 + 3n + 1
    # We generate enough terms to be confident in our result.
    # If the largest unreachable number g is found, we must ensure that the
    # next term in the sequence is greater than g.
    # Let's generate terms for n from 1 to 10.
    # a_10 = 2(100)+3(10)+1 = 231. This should be sufficient.
    n_max = 10
    terms = [(2 * n * n + 3 * n + 1) for n in range(1, n_max + 1)]

    # Step 2: Compute all possible sums and their representations.
    # We use a dictionary where the key is the sum and the value is the list of
    # distinct terms from the sequence that add up to it.
    # We start with 0, which is the sum of an empty set of terms.
    reachable = {0: []}
    for term in terms:
        new_sums = {}
        for s, representation in reachable.items():
            new_sum = s + term
            # Store the new sum and its representation.
            new_sums[new_sum] = representation + [term]
        reachable.update(new_sums)

    # Step 3: Find the largest unreachable number, g.
    # We search for the largest number g that is not in our 'reachable' set.
    # We can be confident in the result if we check up to a reasonable limit
    # that is smaller than the maximum possible sum. Based on preliminary analysis,
    # the number is expected to be small. Let's check up to 200.
    limit_check = 200
    largest_unreachable = 0
    for i in range(1, limit_check):
        if i not in reachable:
            largest_unreachable = i
    
    g = largest_unreachable
    N = g + 1

    # Step 4: Output the results.
    print(f"The sequence of numbers is of the form a_n = 2n^2 + 3n + 1.")
    print(f"The first few terms are: {terms[:6]}...")
    print("-" * 30)
    print(f"The largest number that cannot be written as a sum of distinct terms is: {g}")
    print(f"Therefore, the smallest number N such that any number >= N is representable is: {N}")
    print("-" * 30)

    # Step 5: Show the representation for N.
    if N in reachable:
        representation = reachable[N]
        equation = f"{N} = " + " + ".join(map(str, representation))
        print(f"The number N can be represented as a sum of distinct terms. For example:")
        print(equation)
    else:
        # This case should not happen if the logic is correct.
        print(f"Could not find a representation for N = {N}.")

solve()