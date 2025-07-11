import itertools

# The problem is to find a(n) for n=10.
n = 10

# Initialize a counter for the number of valid ways.
count = 0

# The set of values for y_i is {n+1, ..., 2n}.
y_values = range(n + 1, 2 * n + 1)

# We need to find the number of permutations of y_values that satisfy the condition.
# itertools.permutations generates all possible orderings of y_values.
for y_perm in itertools.permutations(y_values):
    # For each permutation, let's check if it's a valid one.
    # y_perm is a tuple (y_1, y_2, ..., y_n).

    # We use a set to keep track of the generated numbers y_i+i and y_i-i
    # to efficiently check for duplicates.
    generated_values = set()
    is_valid = True

    # We iterate from i = 1 to n.
    for i in range(n):
        # The problem statement uses 1-based indexing for i.
        i_val = i + 1
        
        # y_i is the number from {n+1, ..., 2n} paired with i.
        y_i = y_perm[i]

        # Calculate the two numbers for this pair.
        val_plus = y_i + i_val
        val_minus = y_i - i_val

        # Check if either of the new values already exists in our set of generated values.
        if val_plus in generated_values or val_minus in generated_values:
            # A collision is found, so this permutation is not valid.
            is_valid = False
            break # We can stop checking this permutation and move to the next one.
        
        # If no collision, add the two new unique values to the set.
        generated_values.add(val_plus)
        generated_values.add(val_minus)

    # If the inner loop completed without breaking, is_valid will be true.
    if is_valid:
        # This permutation is valid, so we increment our counter.
        count += 1

# Print the final result in the format a(n) = result.
# This outputs each number in the final equation.
print(f"a({n}) = {count}")