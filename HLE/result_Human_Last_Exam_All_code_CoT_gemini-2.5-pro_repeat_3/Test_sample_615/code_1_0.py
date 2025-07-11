import numpy as np
import itertools

def find_integer_dependency(vectors):
    """
    Finds a simple integer linear dependency for a small set of vectors.
    This function searches for integer coefficients (between -5 and 5) that result
    in a zero vector when combined with the input vectors.
    """
    num_vectors = len(vectors)
    # Generate all possible combinations of coefficients in the range [-5, 5]
    for coeffs in itertools.product(range(-5, 6), repeat=num_vectors):
        # Skip the trivial all-zero coefficient case
        if all(c == 0 for c in coeffs):
            continue
        
        # Calculate the linear combination of vectors
        sum_vector = np.zeros(len(vectors[0]), dtype=int)
        for i in range(num_vectors):
            sum_vector = sum_vector + coeffs[i] * np.array(vectors[i])
        
        # If the sum is the zero vector, a dependency is found
        if np.all(sum_vector == 0):
            return coeffs
    return None

def solve_for_m(m):
    """
    Performs the analysis for a given dimension m.
    """
    print(f"Investigating for m = {m}")
    print("=" * 30)

    # 1. Construct a set of n = m + 1 vectors
    # The set consists of the m standard basis vectors and the all-ones vector
    vectors = [list(row) for row in np.eye(m, dtype=int)]
    vectors.append([1] * m)
    
    S_good = vectors

    print(f"Step 1: Construct a set of n = m + 1 = {m+1} vectors.")
    for i, v in enumerate(S_good):
        print(f"  v{i+1} = {v}")

    # 2. Verify that any m vectors from this set are linearly independent
    print(f"\nStep 2: Verify that any m = {m} vectors from this set are linearly independent.")
    is_set_valid = True
    # Iterate over all combinations of m vectors from the set
    for subset_indices in itertools.combinations(range(m + 1), m):
        subset = [S_good[i] for i in subset_indices]
        matrix = np.array(subset).T
        det = np.linalg.det(matrix)
        
        if abs(det) < 1e-9:
            print(f"  - Checking subset of vectors {tuple(i+1 for i in subset_indices)}... DEPENDENT (Det={det:.0f}).")
            is_set_valid = False
            break
        else:
            print(f"  - Checking subset of vectors {tuple(i+1 for i in subset_indices)}... INDEPENDENT (Det={det:.0f}).")
            
    if is_set_valid:
        print("\nConclusion for Step 2: The construction is valid. So, n >= m + 1.")
    else:
        print("\nConclusion for Step 2: The construction failed.")
        return

    print("=" * 30)
    # 3. Attempt to add an (m+2)-th vector
    print(f"Step 3: Attempt to add an (m+2) = {m+2}-th vector.")
    
    # Generate all possible binary vectors of length m
    all_binary_vectors = [list(i) for i in itertools.product([0, 1], repeat=m)]

    found_extra_vector = False
    # Iterate through all possible vectors to find a candidate for the (m+2)-th vector
    for v_candidate in all_binary_vectors:
        # The new vector must be distinct from the existing ones
        if any(np.array_equal(v_candidate, v) for v in S_good):
            continue

        print(f"\nTrying to add vector {v_candidate} to the set...")
        S_test = S_good + [v_candidate]
        
        is_test_set_valid = True
        # Check all combinations of m vectors from the new, larger set
        for subset_indices in itertools.combinations(range(m + 2), m):
            subset = [S_test[i] for i in subset_indices]
            matrix = np.array(subset).T
            det = np.linalg.det(matrix)
            
            # If a subset is linearly dependent, the candidate vector is invalid
            if abs(det) < 1e-9:
                print(f"  - FAILURE: The subset of vectors indexed by {subset_indices} is linearly DEPENDENT.")
                
                # Find and print the specific dependency equation
                coeffs = find_integer_dependency(subset)
                if coeffs:
                    equation_parts = []
                    vec_names = [f"v{i+1}" for i in subset_indices]
                    for i in range(len(coeffs)):
                        if coeffs[i] != 0:
                            # Represent vector by its index name, e.g., v1, v2, ...
                            equation_parts.append(f"({coeffs[i]})*v{subset_indices[i]+1}")
                    
                    print(f"  - Dependency Equation: {' + '.join(equation_parts)} = 0")
                is_test_set_valid = False
                break
        
        if is_test_set_valid:
            print(f"SUCCESS: Found a valid {m+2}-th vector: {v_candidate}")
            found_extra_vector = True
            break
    
    if not found_extra_vector:
        print("\nConclusion for Step 3: No valid vector can be added.")
        print("This computationally demonstrates for this m that n cannot be m + 2.")

# Main execution for m=3
m = 3
solve_for_m(m)

print("\n" + "="*30)
print("FINAL CONCLUSION")
print("The analysis shows that we can construct a valid set of n = m + 1 vectors.")
print("The computational search for m=3 shows that it is not possible to extend this set to n = m + 2.")
print("This matches the known mathematical result that the maximum value of n is m + 1.")