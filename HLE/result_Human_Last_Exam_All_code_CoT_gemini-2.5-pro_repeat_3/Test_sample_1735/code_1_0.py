import itertools

def get_max_cumulative_impact(seq):
    """Calculates the maximum absolute cumulative sum for a sequence."""
    max_impact = 0
    current_sum = 0
    for x in seq:
        current_sum += x
        if abs(current_sum) > max_impact:
            max_impact = abs(current_sum)
    return max_impact

def count_adjacent_pairs(seq):
    """Counts adjacent pairs of the form {x, -x} in a sequence."""
    count = 0
    for i in range(len(seq) - 1):
        if seq[i] + seq[i+1] == 0:
            count += 1
    return count

def analyze_statement_j():
    """
    Analyzes statement J for a specific counterexample.
    Statement J: For n pairs {x,-x}, at least n-1 pairs are adjacent in any optimal solution.
    """
    A = [1, -1, 2, -2, 3, -3]
    n_pairs = 3
    required_adjacent_pairs = n_pairs - 1

    min_max_impact = float('inf')
    optimal_perms = []

    # Generate all unique permutations of A
    permutations = list(itertools.permutations(A))

    # Find the optimal value by checking all permutations
    for p in permutations:
        impact = get_max_cumulative_impact(p)
        if impact < min_max_impact:
            min_max_impact = impact
            optimal_perms = [p]
        elif impact == min_max_impact:
            optimal_perms.append(p)

    print(f"Analysis for Statement J with input A = {A}")
    print(f"Number of pairs n = {n_pairs}")
    print(f"Optimal maximum cumulative impact found is: {min_max_impact}")
    
    found_counterexample = False
    for p in optimal_perms:
        adj_pairs = count_adjacent_pairs(p)
        if adj_pairs < required_adjacent_pairs:
            print(f"\nFound an optimal permutation that serves as a counterexample:")
            print(f"  Permutation: {p}")
            print(f"  Number of adjacent pairs: {adj_pairs}")
            print(f"  This is less than the required n-1 = {required_adjacent_pairs} pairs.")
            found_counterexample = True
            break  # One counterexample is enough

    if found_counterexample:
        print("\nConclusion: Since an optimal solution exists that violates the condition, Statement J is false.")
    else:
        print("\nConclusion: All optimal solutions satisfy the condition. Statement J holds for this example.")

# Execute the analysis
analyze_statement_j()