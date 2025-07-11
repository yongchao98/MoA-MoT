import math

def solve_task():
    """
    This script calculates the minimal order u_r(n) of the Picard-Fuchs 
    differential equation for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) 
    for n = 3, 4, ..., 12.
    """
    
    print("The minimal order u_r(n) is determined by the symmetry of the potential V(q).")
    print("The degree of the potential is n.")
    print("-" * 40)
    
    results = []
    for n in range(3, 13):
        # Check if n is odd or even to apply the correct formula.
        if n % 2 == 1:
            # For odd n, the potential is not symmetric, so the order is n-1.
            order = n - 1
            print(f"For n = {n} (odd), the formula is u_r(n) = n - 1.")
            print(f"u_r({n}) = {n} - 1 = {order}")
        else:
            # For even n, the potential is an even function, so the order is (n-2)/2.
            order = (n - 2) // 2
            print(f"For n = {n} (even), the formula is u_r(n) = (n - 2) / 2.")
            print(f"u_r({n}) = ({n} - 2) / 2 = {order}")
            
        results.append(order)
        print("-" * 40)

    print("The complete sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    # Using print to output the list directly
    print(results)

solve_task()