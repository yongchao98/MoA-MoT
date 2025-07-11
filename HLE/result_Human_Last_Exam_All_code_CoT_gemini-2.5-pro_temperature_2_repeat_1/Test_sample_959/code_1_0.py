import math

def calculate_sum():
    """
    Calculates the sum described in the problem.

    The logic hinges on the fact that the group G is trivial, which simplifies
    the problem significantly.
    The sum collapses to a single term: the order of the outer automorphism
    group of C_31.
    """

    # The cyclic group is of order p.
    p = 31

    # The set of extensions E has only one element, E = C_p.
    # The sum is the order of Out(E) = Out(C_p).
    # For an abelian group, Out(A) = Aut(A).
    # The order of Aut(C_p) is phi(p).
    # For a prime p, phi(p) = p - 1.
    order_of_out_E = p - 1

    # The sum is over a single term.
    total_sum = order_of_out_E
    
    # The final equation is phi(31) = 31 - 1 = 30.
    # The problem asks to output the numbers in the equation.
    print(f"The calculation is based on Euler's totient function, phi(n).")
    print(f"phi({p}) = {p} - 1 = {total_sum}")
    print(f"The final sum is: {total_sum}")


calculate_sum()
