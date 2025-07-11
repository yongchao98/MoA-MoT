import heapq

def frobenius_number_for_set(numbers):
    """
    Computes the Frobenius number for a set of integers using an algorithm
    based on finding the longest path in a related graph (Dijkstra's).
    The set of numbers must be coprime.
    """
    
    # Sort numbers for convenience, smallest is 'a'
    numbers.sort()
    a = numbers[0]

    # dists[i] will store the smallest number representable by the set
    # that is congruent to i (mod a).
    # We initialize dists for residue 0.
    dists = {0: 0}
    # Priority queue will store (distance, residue) tuples.
    pq = [(0, 0)] 

    while pq:
        dist, r = heapq.heappop(pq)

        # If we've found a shorter path to this residue, skip.
        if dist > dists[r]:
            continue

        # For each number other than 'a', explore new paths.
        for num in numbers[1:]:
            new_r = (r + num) % a
            new_dist = dist + num
            
            # If this is the first path to new_r or a shorter path, update.
            if new_r not in dists or new_dist < dists[new_r]:
                dists[new_r] = new_dist
                heapq.heappush(pq, (new_dist, new_r))
    
    # The Frobenius number is the largest number in these 'longest paths' minus 'a'.
    max_dist = 0
    for r in range(a):
        max_dist = max(max_dist, dists[r])

    return max_dist - a

# The set of numbers determined from the problem analysis
S = [29, 64, 67]

# Calculate the Frobenius number
g = frobenius_number_for_set(S)

# Print the inputs and the final result as requested
print(f"The set of numbers for the Frobenius Coin Problem is {{{S[0]}, {S[1]}, {S[2]}}}.")
print(f"The Frobenius number is {g}.")