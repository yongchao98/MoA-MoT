def count_even_divisors(n):
    """
    Counts the number of even positive divisors of an integer n.
    """
    count = 0
    # Iterate from 1 to n to find all divisors
    for i in range(1, n + 1):
        if n % i == 0:
            # Check if the divisor is even
            if i % 2 == 0:
                count += 1
    return count

# For the generalized quaternion group Q_128, the maximal cyclic subgroup <x> has order 64.
order_of_cyclic_subgroup = 64

# The number of power subgroups coming from even exponents is equal to the
# number of even divisors of the order of the maximal cyclic subgroup.
num_even_power_subgroups = count_even_divisors(order_of_cyclic_subgroup)

# The number of power subgroups coming from odd exponents is 1 (the group Q_128 itself).
num_odd_power_subgroups = 1

# The total number of power subgroups is the sum of the two cases.
total_power_subgroups = num_odd_power_subgroups + num_even_power_subgroups

print("The number of power subgroups in the generalized quaternion group of size 128 can be calculated as follows:")
print("1. From odd powers k, we get one subgroup (the group itself).")
print(f"2. From even powers k, we get a number of subgroups equal to the number of even divisors of {order_of_cyclic_subgroup}.")
print(f"The number of even divisors of {order_of_cyclic_subgroup} is {num_even_power_subgroups}.")
print("\nThe total number of power subgroups is the sum from these two cases:")
print(f"{num_odd_power_subgroups} + {num_even_power_subgroups} = {total_power_subgroups}")