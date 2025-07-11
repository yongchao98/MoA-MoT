import itertools

def calculate_f(a, b, c, d):
    """
    Computes the number of steps for the Ducci sequence to reach (0,0,0,0).
    f(a,b,c,d) includes the initial square, so the count starts at 1.
    """
    current_tuple = (a, b, c, d)
    
    count = 1
    # This process is known to always terminate for integers.
    # We add a safety break for extremely long (or unforeseen cyclic) sequences.
    max_iter = 200
    while current_tuple != (0, 0, 0, 0) and count < max_iter:
        ca, cb, cc, cd = current_tuple
        current_tuple = (abs(ca - cb), abs(cb - cc), abs(cc - cd), abs(cd - ca))
        count += 1
        
    return count

# Step 1: Generate Tribonacci numbers up to the limit 10,000,000
T = [0, 0, 1]
limit = 10000000
while True:
    next_t = T[-1] + T[-2] + T[-3]
    if next_t > limit:
        break
    T.append(next_t)

# The largest n such that T_n+2 <= limit.
# Using 0-based index, T[29] is the largest relevant number.
# Our candidate numbers are T[27], T[28], T[29].
N1 = T[27]  # 2555757
N2 = T[28]  # 4700770
N3 = T[29]  # 8646064

# Step 2: Test all permutations of (0, N1, N2, N3)
base_tuple_elements = [0, N1, N2, N3]
permutations = list(itertools.permutations(base_tuple_elements))

max_len = 0
max_len_tuples = []

for p in permutations:
    length = calculate_f(p[0], p[1], p[2], p[3])
    if length > max_len:
        max_len = length
        max_len_tuples = [p]
    elif length == max_len:
        max_len_tuples.append(p)

# Step 3: Choose the lexicographically smallest tuple from the ones with max length
# This provides a deterministic choice among tuples with the same (minimal) sum.
if not max_len_tuples:
    print("Error: Could not find any tuples.")
else:
    final_tuple = min(max_len_tuples)

    a, b, c, d = final_tuple
    
    # Step 4: Compute the final expression
    result = (a + b - c - d) % 1000

    print(f"The tuple (a, b, c, d) that maximizes the sequence length (M={max_len}) and has the minimum sum is: {final_tuple}")
    print(f"The equation is: a + b - c - d = {a} + {b} - {c} - {d}")
    print(f"The result of (a + b - c - d) mod 1000 is:")
    print(result)
    
    # Final answer in the required format
    # The actual print(result) above gives the number, this is for the specified format.
    # print(f"<<<{result}>>>")