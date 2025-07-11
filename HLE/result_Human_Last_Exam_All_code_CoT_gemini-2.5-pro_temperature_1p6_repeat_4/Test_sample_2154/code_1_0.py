def solve():
    """
    Calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for n from 3 to 12 for the Hamiltonian V(q) = (1/2)(q^2 - q^n).
    """

    print("Calculating the sequence {u_r(3), u_r(4), ..., u_r(12)}:")
    
    results = []
    for n in range(3, 13):
        if n % 2 != 0:  # n is odd
            order = n - 1
            results.append(order)
            print(f"For n = {n} (odd):   u_r({n}) = n - 1 = {n} - 1 = {order}")
        else:  # n is even
            # Using integer division
            order = n // 2
            results.append(order)
            print(f"For n = {n} (even):  u_r({n}) = n / 2 = {n} // 2 = {order}")

    print("\nThe complete sequence is:")
    print(f"{{u_r(3), u_r(4), ..., u_r(12)}} = {{{', '.join(map(str, results))}}}")

solve()