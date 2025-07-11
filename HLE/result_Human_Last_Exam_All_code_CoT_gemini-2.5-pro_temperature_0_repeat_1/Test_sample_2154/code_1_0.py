def solve():
    """
    Calculates the sequence u_r(n) for n from 3 to 12.

    u_r(n) is the minimal order of the Picard-Fuchs differential equation
    for the period of the Hamiltonian with potential V(q) = 1/2 * (q^2 - q^n).

    The order depends on the symmetry of the potential V(q).
    - If n is even, V(q) is an even (symmetric) function, and the order is n/2.
    - If n is odd, V(q) is not a symmetric function, and the order is n-1.
    """
    results = []
    print("Calculating the sequence {u_r(3), u_r(4), ..., u_r(12)}:")
    
    for n in range(3, 13):
        if n % 2 == 0:
            # n is even, potential is symmetric
            order = n // 2
        else:
            # n is odd, potential is not symmetric
            order = n - 1
        
        results.append(order)
        print(f"u_r({n}) = {order}")
        
    print("\nThe complete sequence is:")
    print(results)

solve()