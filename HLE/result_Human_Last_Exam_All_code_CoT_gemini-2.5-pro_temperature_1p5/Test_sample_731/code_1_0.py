import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    For n=2^k, phi(n) = n / 2.
    """
    if n == 1:
        return 1
    # Check if n is a power of 2, which is true for all our orders > 1
    if n > 0 and (n & (n - 1)) == 0:
        return n // 2
    # A general implementation is not strictly necessary here but included for completeness.
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return int(result)

# The semidihedral group is SD_n of order 2^n.
# Order = 512 = 2^9, so n = 9.
n = 9

# N(d): A dictionary to store the number of elements of order d.
N = {}

# Step 1: Count elements from the cyclic subgroup <r> of order 2^(n-1) = 256.
# For each divisor d of 256 (i.e., d=2^k for k=0..8), there are phi(d) elements of order d.
for k in range(n): # k from 0 to 8
    d = 2**k
    N[d] = phi(d)

# Step 2: Add elements of the form r^i * s.
# There are 2^(n-2) = 128 elements of order 2 (where i is even).
# There are 2^(n-2) = 128 elements of order 4 (where i is odd).
num_s_type_half = 2**(n - 2)

# Add 128 elements to the count for order 2
N[2] += num_s_type_half
# Add 128 elements to the count for order 4
N[4] += num_s_type_half

# Step 3: Calculate the number of subgroups for each order and the total.
subgroup_counts_per_order = []
total_subgroups = 0

orders = sorted(N.keys())

print("The number of power subgroups is the sum of the number of cyclic subgroups of each possible order.")
print("Number of subgroups of order d = (Number of elements of order d) / phi(d)\n")

for d in orders:
    elements_of_order_d = N[d]
    phi_d = phi(d)
    subgroups_of_order_d = elements_of_order_d // phi_d
    
    subgroup_counts_per_order.append(subgroups_of_order_d)
    
    print(f"Order d={d}:")
    print(f"  Number of elements = {elements_of_order_d}")
    print(f"  phi({d}) = {phi_d}")
    print(f"  Number of subgroups = {elements_of_order_d} / {phi_d} = {subgroups_of_order_d}\n")
    
total_subgroups = sum(subgroup_counts_per_order)

# Step 4: Print the final equation and result, as requested.
equation_str = " + ".join(map(str, subgroup_counts_per_order))
print(f"The final calculation is the sum of the number of subgroups for each order:")
print(f"{equation_str} = {total_subgroups}")

print("\nTherefore, the total number of power subgroups is:")
print(total_subgroups)