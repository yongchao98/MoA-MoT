def solve_exponents():
    """
    This function calculates the value of sum(k * n_k) for k=1 to 6
    based on the derived exponents for the magnetic field noise problem.
    """

    # Exponents from the low-frequency (lim omega -> 0) scaling law:
    # S_B(omega) proportional to sigma^(n1) * T^(n2) * z^(n3)
    n1 = 1
    n2 = 1
    n3 = -2

    # Exponents from the frequency spectrum scaling S_B(omega) ~ omega^n
    # in different regimes.
    n4 = 0    # for omega << 1/(sigma*z*t)
    n5 = -2   # for 1/(sigma*z*t) << omega << 1/(sigma*t^2)
    n6 = -1   # for omega >> 1/(sigma*t^2)
    
    exponents = [n1, n2, n3, n4, n5, n6]
    
    # We want to calculate the sum of k * n_k for k from 1 to 6.
    # The list 'exponents' holds n_1, n_2, ..., n_6 at indices 0, 1, ..., 5.
    # So we need to calculate sum_{i=0 to 5} (i+1) * exponents[i]
    
    total_sum = 0
    equation_parts = []
    
    for k in range(1, 7):
        n_k = exponents[k-1]
        term = k * n_k
        total_sum += term
        # Store each part of the equation for printing
        equation_parts.append(f"({k} * {n_k})")
        
    equation_str = " + ".join(equation_parts)
    print(f"Calculating the sum of k * n_k for k=1 to 6:")
    print(f"The exponents are: n1={n1}, n2={n2}, n3={n3}, n4={n4}, n5={n5}, n6={n6}")
    print(f"The calculation is: {equation_str}")
    print(f"Result of each term: {[k * exponents[k-1] for k in range(1,7)]}")
    print(f"The final sum is: {total_sum}")

solve_exponents()