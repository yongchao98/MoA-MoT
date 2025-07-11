import collections

def get_divisors(n):
    """Returns a sorted list of all divisors of n."""
    divs = set()
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def find_subset_sum(numbers, target):
    """
    Checks if a subset of 'numbers' can sum to 'target' and returns the subset.
    Uses a dynamic programming approach with backtracking to find the actual numbers.
    """
    dp = {0: []}
    for num in numbers:
        # Create new sums by adding the current number to existing sums
        new_sums = {}
        for s, subset in dp.items():
            if s + num <= target:
                new_sums[s + num] = subset + [num]
        
        # Add the new sums to the dp table
        for s, subset in new_sums.items():
            # Only update if the key is new, to keep the first found subset
            if s not in dp:
                dp[s] = subset

        if target in dp:
            return dp[target]
            
    return None if target not in dp else dp[target]

def is_zumkeller(n):
    """
    Checks if a number n is a Zumkeller number.
    Returns a tuple (is_zumkeller, details_dictionary).
    """
    divisors = get_divisors(n)
    sigma = sum(divisors)

    # Sum of divisors must be even
    if sigma % 2 != 0:
        return False, {}

    target_partition_sum = sigma // 2

    # n itself must be a divisor. If n > target_partition_sum, it can't be
    # in any partition, so it's not a Zumkeller number.
    if n > target_partition_sum:
        return False, {}

    # We need to find a subset of the other divisors that sums to the remainder
    proper_divisors = [d for d in divisors if d != n]
    remainder_needed = target_partition_sum - n
    
    subset = find_subset_sum(proper_divisors, remainder_needed)

    if subset is not None:
        details = {
            "n": n,
            "sigma": sigma,
            "target_sum": target_partition_sum,
            "remainder": remainder_needed,
            "subset": subset
        }
        return True, details
    
    return False, {}

def solve():
    """Finds the next number in the sequence."""
    last_num_in_sequence = 2352
    current_num = last_num_in_sequence + 1

    while True:
        is_zum, details = is_zumkeller(current_num)
        if is_zum:
            # Found the next number
            print(f"The sequence is composed of Zumkeller numbers.")
            print(f"The last number given is {last_num_in_sequence}.")
            print(f"The next Zumkeller number is {details['n']}.\n")
            print("Here is the proof:")
            print(f"1. The sum of all divisors of {details['n']} is {details['sigma']}.")
            print(f"2. The divisors must be split into two sets, each summing to {details['sigma']} / 2 = {details['target_sum']}.")
            print(f"3. One partition includes {details['n']}. The other divisors in this partition must sum to {details['target_sum']} - {details['n']} = {details['remainder']}.")
            
            # Formatting the final equation with each number
            subset_str = ' + '.join(map(str, details['subset']))
            print(f"4. We found a subset of other divisors that fits: {subset_str} = {details['remainder']}.")
            print(f"\nTherefore, the final answer for the next number in the sequence is {details['n']}.")
            break
        current_num += 1

solve()
<<<2360>>>