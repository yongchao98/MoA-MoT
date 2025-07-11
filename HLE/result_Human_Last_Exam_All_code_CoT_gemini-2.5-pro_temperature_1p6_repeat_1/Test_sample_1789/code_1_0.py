import math

# The upper bound of the interval.
N = 823564528378596

# The number of antichains is 20.
k = 20

# We will store each term of the sum to display the full equation.
terms = []
total_size = 0

# The simplified formula for the total size is:
# Sum_{v=0 to k-1} (number of odd integers <= N / 2^v)
# The number of odd integers <= L is (L + 1) // 2.
for v in range(k):
    # Calculate the limit for the odd numbers d, where d * 2^v <= N.
    limit = N // (2**v)
    
    # Calculate the number of odd integers up to this limit.
    num_odd = (limit + 1) // 2
    
    terms.append(num_odd)
    total_size += num_odd

# Format the final output to show the full equation as requested.
equation_str = " + ".join(map(str, terms))
print(f"{equation_str} = {total_size}")