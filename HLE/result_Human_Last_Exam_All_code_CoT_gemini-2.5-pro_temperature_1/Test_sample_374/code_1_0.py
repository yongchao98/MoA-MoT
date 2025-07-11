import math

def get_gl_order(n, q):
    """Computes the order of the general linear group GL(n, q)."""
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def get_odd_part(num):
    """Returns the largest odd divisor of a number."""
    while num % 2 == 0:
        num //= 2
    return num

def solve():
    """
    Solves the problem by following the outlined plan.
    """
    n = 4
    q = 2
    
    # Step 1 & 2: E is an odd-order subgroup of GL(4, F_2)
    print(f"The inertial quotient E is an odd-order subgroup of Aut(D) = GL({n}, F_{q}).")
    gl_order = get_gl_order(n, q)
    print(f"The order of GL({n}, F_{q}) is (2^{n}-1)(2^{n}-2)(2^{n}-4)(2^{n}-8) = 15 * 14 * 12 * 8 = {gl_order}.")
    
    odd_order_bound = get_odd_part(gl_order)
    print(f"The maximum possible order for an odd-order subgroup must divide the odd part of |GL({n}, F_{q})|, which is {odd_order_bound}.\n")
    
    # Step 3, 4, 5: E must be cyclic and act freely.
    print("A crucial constraint from block theory is that E must act freely on the non-zero elements of D.")
    print("This implies that no non-identity element of E can have 1 as an eigenvalue.")
    print("A group with such a linear action must be cyclic. So, E is a cyclic group.")
    print("We need to find the maximum odd order of an element g in GL(4, F_2) whose minimal polynomial is not divisible by (x+1).\n")
    
    # Step 6: Analyze element orders based on characteristic polynomials.
    print("Analyzing the orders of elements based on their characteristic polynomials (degree 4):")
    
    # Case 1: The characteristic polynomial is irreducible of degree 4.
    # These polynomials do not have (x+1) as a factor.
    print("  Case 1: The characteristic polynomial is irreducible.")
    # The order of an element whose characteristic polynomial is an irreducible poly of degree n over F_q
    # is the order of the multiplicative group of the field F_{q^n}, which is q^n-1, or a divisor of it.
    # For x^4+x+1, the order is 2^4-1 = 15.
    order1 = 15
    print(f"    - For a polynomial like x^4+x+1, the element order is {order1}.")
    # For x^4+x^3+x^2+x+1 = (x^5-1)/(x-1), the order is 5.
    order2 = 5
    print(f"    - For the polynomial x^4+x^3+x^2+x+1, the element order is {order2}.")
    
    # Case 2: The characteristic polynomial is reducible.
    # It must be a product of irreducibles, none of which is (x+1).
    # The only other way to partition 4 is 2+2.
    # The only irreducible polynomial of degree 2 over F_2 is x^2+x+1.
    print("\n  Case 2: The characteristic polynomial is reducible.")
    print("    The only possible factorization without (x+1) as a factor is (x^2+x+1)^2.")
    # The minimal polynomial is x^2+x+1. The order of the corresponding element is 3 (2^2-1).
    order3 = 3
    print(f"    - This corresponds to an element of order {order3}.")
    
    max_order = max(order1, order2, order3)
    print(f"\nThe possible orders for such an element g are divisors of {order1}, {order2}, or {order3}.")
    print(f"The maximum of these possible orders is {max_order}.")
    print("\nA cyclic group of order 15 exists in GL(4, F_2) and satisfies the free-action condition.")
    print(f"\nTherefore, the highest possible order that E can have is {max_order}.")

solve()