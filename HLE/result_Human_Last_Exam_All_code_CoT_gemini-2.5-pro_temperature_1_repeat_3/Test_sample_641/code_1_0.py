def solve():
    """
    Calculates the number of involutions in the group PSU(4, 997).
    """
    q = 997

    # Calculate the size of the first class of involutions (N1)
    # These arise from elements g in SU(4,q) with g^2 = I.
    # Formula: N1 = q^4 * (q^2 - q + 1) * (q^2 + 1)
    q2 = q * q
    q4 = q2 * q2
    
    n1_factor1 = q2 - q + 1
    n1_factor2 = q2 + 1
    N1 = q4 * n1_factor1 * n1_factor2

    # Calculate the size of the second class of involutions (N2)
    # These arise from elements g in SU(4,q) with g^2 = -I.
    # Formula: N2 = q^4 * (q + 1) * (q^3 + 1)
    q3 = q2 * q
    
    n2_factor1 = q + 1
    n2_factor2 = q3 + 1
    N2 = q4 * n2_factor1 * n2_factor2
    
    # The total number of involutions is the sum of the sizes of these two classes.
    total_involutions = N1 + N2
    
    # Print the full equation with all numbers
    print(f"The number of involutions in PSU(4, 997) is the sum of the sizes of two classes, N1 and N2.")
    print(f"N1 = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1) = {q4} * {n1_factor1} * {n1_factor2} = {N1}")
    print(f"N2 = {q}^4 * ({q} + 1) * ({q}^3 + 1) = {q4} * {n2_factor1} * {n2_factor2} = {N2}")
    print(f"Total = N1 + N2 = {N1} + {N2} = {total_involutions}")

solve()