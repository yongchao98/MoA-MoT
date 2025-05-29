import math

# List of pairs
pairs = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), 
         (2, 3), (2, 4), (2, 5), (2, 6), 
         (3, 4), (3, 5), (3, 6), 
         (4, 5), (4, 6), 
         (5, 6)]

# Count coprime pairs
coprime_count = sum(1 for a, b in pairs if math.gcd(a, b) == 1)

# Total number of pairs
total_pairs = len(pairs)

# Probability
probability = coprime_count / total_pairs

print(f"{coprime_count}/{total_pairs}")