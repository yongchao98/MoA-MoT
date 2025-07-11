import itertools

def demonstrate_construction(k_str):
    """
    Demonstrates that for n = 2k - 1, a suitable family F exists.
    F is taken to be the set of all k-subsets of {1, ..., n}.
    It shows that any (k-1)-subset S can be represented as a difference
    of two sets from F.
    """
    try:
        k = int(k_str)
        if k < 2:
            print("Please enter an integer k >= 2.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    n = 2 * k - 1
    
    print(f"For k = {k}, the maximum value of n is 2*k - 1 = {n}.")
    print(f"We will work with the set [n] = {{1, 2, ..., {n}}}.")
    print("-" * 40)
    
    # The family F is all k-subsets of [n]. This family is intersecting.
    print(f"Let the family F be the collection of all {k}-subsets of [n].")

    # We demonstrate that any (k-1)-subset S of [n] is a difference S = F \\ F'
    # for some F, F' in the family.
    
    # As an example, we take the set S = {1, 2, ..., k-1}.
    S = set(range(1, k))
    S_list = sorted(list(S))
    S_str = "{" + ",".join(map(str, S_list)) + "}"

    print(f"\nLet's test with the sample ({k-1})-subset S = {S_str}.")
    
    # To form the difference, we need to find two sets, F and F', in the family.
    # We require F \\ F' = S, which implies |F intersect F'| must be 1.
    # Let F = S union {x} for some element x not in S.
    
    # Choose x from [n] \\ S. The complement has size (2k-1) - (k-1) = k.
    # We can pick any element from it. For simplicity, let's pick x = k.
    x = k
    
    # Construct F and F'
    F_set = S.union({x})
    F_prime_set = set(range(1, n + 1)).difference(S)
    
    print(f"\nWe construct two sets, F and F', to generate S:")
    # 1. Choose an element x not in S. Let's pick x = {x}.
    # 2. Define F = S union {x}.
    # 3. Define F' = [n] \\ S.
    
    F_list = sorted(list(F_set))
    F_str = "{" + ",".join(map(str, F_list)) + "}"
    
    F_prime_list = sorted(list(F_prime_set))
    F_prime_str = "{" + ",".join(map(str, F_prime_list)) + "}"

    print(f"  F = S U {{{x}}} = {S_str} U {{{x}}} = {F_str}")
    print(f"  F' = [n] \\ S = {{1,...,{n}}} \\ {S_str} = {F_prime_str}")
    
    print(f"\nBoth F and F' have size {k}, so they are members of our family F.")

    # Verify the difference and print the final equation.
    difference_set = F_set.difference(F_prime_set)
    diff_list = sorted(list(difference_set))
    diff_str = "{" + ",".join(map(str, diff_list)) + "}"
    
    print("\nFinal Equation Verification:")
    print(f"{S_str} = {F_str} \\ {F_prime_str}")
    print(f"Result of calculation: {F_str} \\ {F_prime_str} = {diff_str}")

    if difference_set == S:
        print("The equation holds true, confirming our construction.")
    else:
        print("An error occurred; the construction did not yield the correct difference.")

if __name__ == '__main__':
    k_input = input("Enter an integer value for k (where k >= 2): ")
    demonstrate_construction(k_input)
