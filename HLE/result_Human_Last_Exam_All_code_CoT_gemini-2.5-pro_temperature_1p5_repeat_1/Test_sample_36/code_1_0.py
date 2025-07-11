import math

def get_divisors(n):
    """Helper function to find all divisors of an integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return divs

def solve_group_representation_problem():
    """
    Calculates the number of unique eigenvalue sets for representations
    of Abelian groups of order 18.
    """
    order = 18
    print(f"This script finds the number of unique eigenvalue sets for representations of Abelian groups of order {order}.\n")
    
    print("Step 1: Identify the Abelian groups of order 18.")
    print("The two non-isomorphic Abelian groups of order 18 are:")
    print("  G1 = Z_18 (the cyclic group)")
    print("  G2 = Z_2 x Z_3 x Z_3\n")

    print("Step 2: Determine the possible eigenvalue sets for each group.")
    print("For any representation rho, the eigenvalue set S(rho) is a set of m-th roots of unity.")
    print("The possible values for m are the divisors of the exponent of the group.\n")

    # Exponent and divisors for G1
    exp_g1 = 18
    divisors_g1 = get_divisors(exp_g1)
    print(f"For G1 = Z_18, the exponent is {exp_g1}.")
    print(f"The set of possible values for m is D({exp_g1}) = {sorted(list(divisors_g1))}.\n")

    # Exponent and divisors for G2
    # In Python 3.9+, math.lcm can take multiple arguments.
    # For compatibility, we can calculate it iteratively if needed.
    exp_g2 = math.lcm(2, 3, 3)
    divisors_g2 = get_divisors(exp_g2)
    print(f"For G2 = Z_2 x Z_3 x Z_3, the exponent is lcm(2, 3, 3) = {exp_g2}.")
    print(f"The set of possible values for m is D({exp_g2}) = {sorted(list(divisors_g2))}.\n")
    
    print("Step 3: Calculate the total number of unique sets.")
    print("The total number of unique sets is the size of the union of the sets of possible m-values from each group.\n")
    
    all_possible_m = divisors_g1.union(divisors_g2)
    
    print("Calculation:")
    print(f"  | D({exp_g1}) U D({exp_g2}) |")
    print(f"= | {sorted(list(divisors_g1))} U {sorted(list(divisors_g2))} |")
    print(f"= | {sorted(list(all_possible_m))} |")
    
    final_count = len(all_possible_m)
    print(f"= {final_count}\n")
    
    print(f"Thus, there are {final_count} unique sets.")

# Run the solver
solve_group_representation_problem()

# The final answer in the requested format
final_answer = len(get_divisors(18).union(get_divisors(6)))
print(f"<<<{final_answer}>>>")