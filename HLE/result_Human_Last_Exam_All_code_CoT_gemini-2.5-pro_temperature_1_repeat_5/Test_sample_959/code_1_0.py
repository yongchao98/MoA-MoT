def main():
    """
    This script calculates the required sum based on the properties of the group G and its extensions.
    The reasoning is as follows:
    1. The group G is shown to be the trivial group {1}.
    2. Consequently, there is only one central extension E of G by C_31, and E is isomorphic to C_31.
    3. The sum is therefore just the order of the outer automorphism group of C_31.
    4. For an abelian group like C_31, this is the order of its automorphism group, which is phi(31).
    """

    # The order of the cyclic group C is 31.
    n = 31
    
    # The sum is over a single extension E ~= C_n.
    # The value to be calculated is o(E) = |Out(C_n)|.
    # For an abelian group, Out(G) = Aut(G).
    # The order of Aut(C_n) is given by Euler's totient function, phi(n).
    
    # Since n=31 is a prime number, phi(n) = n - 1.
    order_out_E = n - 1

    # The final sum is simply this value.
    total_sum = order_out_E
    
    # The problem asks to output each number in the final equation.
    # The final equation is: sum = o(E) = phi(31) = 31 - 1 = 30.
    print(f"The calculation is: {total_sum} = {n} - 1")
    print(f"The final sum is: {total_sum}")

if __name__ == "__main__":
    main()