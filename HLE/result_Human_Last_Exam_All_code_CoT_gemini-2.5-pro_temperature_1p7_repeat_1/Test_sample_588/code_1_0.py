def solve():
    """
    Solves the problem and prints the result as specified.
    The solution is derived mathematically and this code serves to verify
    and format the output according to the instructions.
    """
    # The lexicographically least tuple derived from the analysis.
    result_tuple = (0, 0, 0, 0, 0, 2)
    
    l = len(result_tuple) // 2
    pairs = []
    for i in range(l):
        pairs.append((result_tuple[2*i], result_tuple[2*i+1]))
        
    a1, b1 = pairs[0]
    a2, b2 = pairs[1]
    a3, b3 = pairs[2]

    # The condition for the connected sum N to be full is that its Euler
    # characteristic is zero.
    # chi(N) = chi(M1) + chi(M2) + chi(M3) - 2*(l-1) = 0
    
    # Calculate chi for each component manifold M(a_i, b_i)
    # chi(M(a,b)) = 4 * (1-a) * (1-b)
    chi1 = 4 * (1 - a1) * (1 - b1)
    chi2 = 4 * (1 - a2) * (1 - b2)
    chi3 = 4 * (1 - a3) * (1 - b3)
    
    # The final equation demonstrates that the condition is met.
    # The sum of the Euler characteristics of the component manifolds
    # must equal 2*(l-1) for the connected sum to be full.
    total_chi = chi1 + chi2 + chi3
    l_minus_1_doubled = 2 * (l - 1)
    
    # Each manifold M(a_i,b_i) is not full, as a_i, b_i are not 1, so chi_i != 0.
    # The connected sum is full, as the following equation holds.
    # I am outputting each number in the final equation.
    print(f"{chi1} + {chi2} + {chi3} - {l_minus_1_doubled} = {total_chi - l_minus_1_doubled}")

    # Final answer formatted as requested.
    flat_tuple_str = str(result_tuple).replace(" ", "")
    print(flat_tuple_str)

solve()
<<< (0,0,0,0,0,2) >>>