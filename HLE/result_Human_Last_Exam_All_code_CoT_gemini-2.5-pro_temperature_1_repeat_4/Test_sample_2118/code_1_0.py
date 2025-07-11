import math

# Set the limit for the exponent n
limit = 100

# --- Case 1: n = 2^p, where p is odd ---
# We count the number of n <= 100 that are powers of 2 with an odd exponent.
cat1_count = 0
cat1_p_values = []
p = 1
while True:
    n = 2**p
    if n > limit:
        break
    cat1_count += 1
    cat1_p_values.append(p)
    p += 2

# --- Case 2: n = m * 2^p, where m > 1 is odd and p is even ---
# We count the number of n <= 100 for each valid even p.
cat2_counts_by_p = {}
p = 0
while True:
    power_of_2 = 2**p
    # The smallest n in this category for a given p is 3 * 2^p.
    if 3 * power_of_2 > limit:
        break
    
    # We need to find the number of odd integers m >= 3 such that m * 2^p <= limit.
    # This is equivalent to m <= limit / 2^p.
    max_m = limit // power_of_2
    
    # If the maximum possible m is less than 3, there are no terms for this p.
    if max_m >= 3:
        # We count the odd numbers in the range [3, max_m].
        # If max_m is even, the largest odd number is max_m - 1.
        if max_m % 2 == 0:
            max_m -= 1
        current_p_count = (max_m - 3) // 2 + 1
        cat2_counts_by_p[p] = current_p_count

    p += 2

# --- Print the results and the final equation ---
print("The number of non-zero terms is calculated by considering two cases based on the decomposition n = m * 2^p (m odd):")

print(f"\nCase 1: n = 2^p (m=1). The coefficient is non-zero if p is odd.")
print(f"For n <= {limit}, the valid odd exponents are p = {', '.join(map(str, cat1_p_values))}, corresponding to n = 2, 8, 32.")
print(f"This gives {cat1_count} terms.")

print(f"\nCase 2: n = m * 2^p (m > 1 odd). The coefficient is non-zero if p is even.")
cat2_total_count = sum(cat2_counts_by_p.values())
cat2_counts_list = list(cat2_counts_by_p.values())

for p_val, count in cat2_counts_by_p.items():
    print(f" - For p = {p_val}, there are {count} valid odd values for m.")

print(f"The total number of terms from Case 2 is the sum: { ' + '.join(map(str, cat2_counts_list)) } = {cat2_total_count}.")

total_count = cat1_count + cat2_total_count
print(f"\nThe total number of non-zero terms is the sum of the counts from both cases.")
print(f"Final Equation: {cat1_count} + {cat2_total_count} = {total_count}")