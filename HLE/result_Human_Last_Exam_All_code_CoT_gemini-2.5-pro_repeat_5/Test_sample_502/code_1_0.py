import math

def find_max_dimension():
    """
    This script determines the largest possible dimension for R/I by investigating
    the existence of a suitable complex reflection group, based on the principles
    of invariant theory.
    """
    n = 10
    group_order = 10000

    print("Step 1: Analyze irreducible complex reflection groups of the G(m, p, n) family.")
    print(f"The group G acts on a space of dimension n = {n}.")
    print(f"The required order of the group is |G| = {group_order}.")
    
    # The order of a group from the infinite family G(m, p, n) is given by the formula:
    # |G(m, p, n)| = (m^n * n!) / p, where p is a positive integer that divides m.
    n_factorial = math.factorial(n)
    print(f"The formula for the order of G(m, p, {n}) is (m^{n} * n!) / p.")
    print(f"For n=10, this is (m^10 * {n_factorial}) / p = {group_order}.")
    
    # We can rearrange this to solve for p: p = m^10 * (n! / |G|)
    factor = n_factorial / group_order
    print(f"Rearranging the formula gives: p = m^10 * {factor}.")
    print("For m and p to be integers such that p divides m, this equation has no solution.")
    print("This implies that no irreducible complex reflection group from this family meets the criteria.\n")

    print("Step 2: Analyze reducible complex reflection groups.")
    # A reducible reflection group can be formed as a direct product of reflection groups
    # acting on smaller dimensional subspaces. The simplest such construction is a product
    # of cyclic groups, C_{m_i}, each acting on a 1D space C^1.
    # The group G = C_{m_1} x ... x C_{m_10} acts on C^10 and is a reflection group.
    # The order of this group is |G| = m_1 * m_2 * ... * m_10.
    print(f"We check if a reducible reflection group on C^{n} of order {group_order} can be constructed.")
    print(f"Consider a group G formed by a product of {n} cyclic groups, G = C_m_1 x ... x C_m_{n}.")
    print(f"The order of G is the product of the orders of the cyclic groups: |G| = m_1 * ... * m_{n}.")
    print(f"We need to find positive integers m_i whose product is {group_order}.")
    
    # We can provide a simple example of such integers.
    m_factors = [1] * n
    m_factors[0] = group_order
    
    product = 1
    for factor_val in m_factors:
        product *= factor_val
        
    print(f"An example solution is m = {m_factors}.")
    
    # Building the string for the equation
    equation_str_parts = [str(f) for f in m_factors]
    print(f"Checking the product: {' * '.join(equation_str_parts)} = {product}.")

    if product == group_order:
        print(f"\nA valid set of integers exists. For instance, the cyclic group C_{group_order} acting on")
        print("the first coordinate of C^10 is a complex reflection group of the required order.")
        print("This confirms that a complex reflection group of order 10000 acting on C^10 exists.\n")
    else:
        # This case is not expected to be reached.
        print("No simple solution found for the product of integers.")

    print("Step 3: Conclusion based on Invariant Theory.")
    # The theorem states dim(R/I) <= |G|, with equality for reflection groups.
    # Since we have demonstrated that such a group exists, the maximum dimension is |G|.
    max_dim = group_order
    print("By the Chevalley-Shephard-Todd theorem, the dimension of the coinvariant algebra R/I")
    print("is equal to the order of the group G if, and only if, G is a complex reflection group.")
    print("Since such a group exists, the largest possible dimension for R/I is the order of the group.")
    
    # Final equation output, as requested.
    print("\nThe final equation for the largest possible dimension is:")
    final_equation_var = "max_dim(R/I)"
    final_equation_val1 = "|G|"
    final_equation_val2 = max_dim
    print(f"{final_equation_var} = {final_equation_val1}")
    print(f"{final_equation_var} = {final_equation_val2}")

find_max_dimension()