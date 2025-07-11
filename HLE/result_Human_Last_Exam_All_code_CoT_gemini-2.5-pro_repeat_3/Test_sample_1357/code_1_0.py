def count_stable_reduction_types():
    """
    This script calculates the number of types of stable reductions of a genus 4 curve
    whose Jacobian has good reduction.
    """

    # The genus of the curve in question.
    g = 4

    print("This problem can be solved by applying a key theorem from arithmetic geometry.")
    print("The theorem states: A curve's Jacobian has good reduction if and only if the curve's stable reduction is irreducible (i.e., it has only one component).")
    print("\nTherefore, the problem is reduced to counting the number of types of irreducible stable curves of arithmetic genus 4.")
    print("-" * 50)

    print(f"An irreducible stable curve of arithmetic genus g = {g} is classified by its number of nodes, denoted by 'δ'.")
    print("The number of nodes determines the geometric genus of the curve's normalization via the formula: g_normalization = g - δ.")
    
    print("\nSince the geometric genus of any curve must be non-negative (>= 0), we get the following constraint:")
    print(f"  g - δ >= 0   =>   {g} - δ >= 0   =>   δ <= {g}")

    print(f"\nThis means the number of nodes 'δ' can be any integer from 0 to {g}.")
    possible_delta_values = list(range(g + 1))
    print(f"The possible values for δ are: {possible_delta_values}.")
    print("-" * 50)
    
    print("Each possible value for δ corresponds to a unique type of stable reduction:")
    print(f" δ=0: A smooth curve of genus {g-0}")
    print(f" δ=1: An irreducible curve with 1 node (normalization has genus {g-1})")
    print(f" δ=2: An irreducible curve with 2 nodes (normalization has genus {g-2})")
    print(f" δ=3: An irreducible curve with 3 nodes (normalization has genus {g-3})")
    print(f" δ=4: An irreducible curve with 4 nodes (normalization has genus {g-4})")
    print("-" * 50)

    # The final calculation is the total count of these possibilities.
    num_types = len(possible_delta_values)
    
    print("The final calculation is the count of possible types, which is based on the genus.")
    print(f"Number of Types = (genus) + 1")
    # As requested, printing the numbers in the final equation
    print(f"Final Equation: {g} + 1 = {num_types}")

count_stable_reduction_types()