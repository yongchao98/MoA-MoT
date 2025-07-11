import math

def phi(n):
    """
    Calculates Euler's totient function for n.
    For this problem, we only need it for powers of 2.
    phi(1) = 1
    phi(2^k) = 2^(k-1) for k > 0
    """
    if n == 1:
        return 1
    # A simple check for powers of 2
    if n > 0 and (n & (n - 1)) == 0:
        return n // 2
    # Fallback for general integers (not needed for this specific problem)
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

def count_power_subgroups_q128():
    """
    Calculates and explains the number of power subgroups in the
    generalized quaternion group of size 128 (Q_128).
    """
    print("Finding the number of power subgroups (cyclic subgroups) in Q_128.")
    print("The strategy is to count the number of elements of each possible order and then determine the number of subgroups for each order.")
    print("\nStep 1: Count the number of elements of each order.")
    
    # Orders of elements in the cyclic subgroup <x> of order 64
    elements_from_x = {
        1: phi(1),   # Order 1:  1 element (x^0)
        2: phi(2),   # Order 2:  1 element (x^32)
        4: phi(4),   # Order 4:  2 elements
        8: phi(8),   # Order 8:  4 elements
        16: phi(16), # Order 16: 8 elements
        32: phi(32), # Order 32: 16 elements
        64: phi(64)  # Order 64: 32 elements
    }
    
    # The other 64 elements are of the form x^k*y. All of them have order 4.
    num_xky_elements = 64
    
    # Total elements of order 4
    total_elements_order_4 = elements_from_x[4] + num_xky_elements

    print(" - Elements from the cyclic part <x>:")
    for order, count in elements_from_x.items():
        print(f"   - Order {order}: {count} elements")

    print(f" - The other {num_xky_elements} elements (of the form x^k*y) all have order 4.")
    
    # Compiling the counts for the entire group
    elements_by_order = {
        1: elements_from_x[1],
        2: elements_from_x[2],
        4: total_elements_order_4,
        8: elements_from_x[8],
        16: elements_from_x[16],
        32: elements_from_x[32],
        64: elements_from_x[64]
    }
    
    print("\nStep 2: Calculate the number of cyclic subgroups for each order.")
    print("This is (Number of elements of order d) / phi(d).")
    
    subgroups_by_order = {}
    for order, num_elements in elements_by_order.items():
        num_subgroups = num_elements // phi(order)
        subgroups_by_order[order] = num_subgroups
        print(f" - Order {order}: {num_elements} elements / phi({order})={phi(order)} -> {num_subgroups} subgroup(s)")
        
    print("\nStep 3: Sum the results to find the total.")
    
    total_subgroups = sum(subgroups_by_order.values())
    
    equation_parts = [str(count) for order, count in sorted(subgroups_by_order.items())]
    equation_str = " + ".join(equation_parts)
    
    print(f"The final calculation is:")
    print(f"{equation_str} = {total_subgroups}")

if __name__ == '__main__':
    count_power_subgroups_q128()
    # The total number is 1+1+33+1+1+1+1 = 39
    print("\n<<<39>>>")
