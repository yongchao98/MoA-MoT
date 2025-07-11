import math

def check_filled_groups_of_order_2q_m(q: int, m: int):
    """
    This function explains why all nonabelian groups of order 2*q**m 
    are filled, for a given odd prime q and natural number m. It does
    so by walking through a proof by contradiction.

    Args:
        q: An odd prime number.
        m: A natural number (positive integer).
    """
    if not isinstance(q, int) or q <= 1 or (q > 2 and q % 2 == 0):
        print(f"Error: q must be an odd prime, but got {q}.")
        return
    if not isinstance(m, int) or m < 1:
        print(f"Error: m must be a natural number (integer >= 1), but got {m}.")
        return

    order = 2 * (q**m)
    q_subgroup_order = q**m

    print(f"Analysis for groups of order 2 * {q}^{m} = {order}")
    print("="*50)
    print("Based on a known theorem, a group G of this order is 'filled' if its unique normal subgroup Q of order q^m is not a union of two of its own conjugacy classes.")
    print("We will prove that Q can never be a union of two conjugacy classes.\n")

    print("--- The Proof by Contradiction ---")
    print(f"1. Let's assume the subgroup Q (of order q^m = {q_subgroup_order}) IS a union of two conjugacy classes, C1 and C2.")
    print("2. The identity element {1} is always a conjugacy class of size 1. So, one class, say C1, must be {{1}}.")
    print("3. This implies the other class, C2, must contain all other elements of Q.")
    
    size_of_c2 = q_subgroup_order - 1
    print(f"4. The size of C2 would be |Q| - 1 = {q_subgroup_order} - 1 = {size_of_c2}.")
    
    print("\n5. A core group theory fact: The size of a conjugacy class must divide the order of the group.")
    print(f"   This means |C2| must divide |Q|. So, {size_of_c2} must divide {q_subgroup_order}.")

    # This is the "final equation" part requested by the user prompt
    print("\n6. To check if this is possible, we examine the final equation from the divisibility rule:")
    print(f"   Equation: ({q_subgroup_order} - 1) | {q_subgroup_order}")
    print(f"   Let's check the greatest common divisor (GCD) of these two numbers:")
    
    # Python code outputting the numbers in the equation
    print(f"   GCD({q_subgroup_order}, {size_of_c2})")

    gcd_result = math.gcd(q_subgroup_order, size_of_c2)
    # Python code outputting the result
    print(f"   The result is {gcd_result}.")
    
    print(f"\n7. Since the GCD is 1, the only way for {size_of_c2} to divide {q_subgroup_order} is if {size_of_c2} = 1.")
    print(f"   If {size_of_c2} = 1, then {q_subgroup_order} = 2.")
    print(f"8. This leads to the equation q^m = {q}^{m} = 2.")
    print(f"   This is a CONTRADICTION, because q is an odd prime ({q}) and m is a natural number ({m}).\n")

    print("--- Conclusion ---")
    print("The initial assumption was false. A q-group Q can never be a union of two conjugacy classes.")
    print("Therefore, by the theorem, EVERY group of order 2*q^m is a filled group.")
    print("\nThe answer to the question 'What are the nonabelian filled groups of order 2q^m?' is:")
    print("ALL nonabelian groups of order 2q^m are filled.")

# You can test with any odd prime q and natural number m.
# Example: q=3, m=2 (groups of order 18)
check_filled_groups_of_order_2q_m(q=3, m=2)