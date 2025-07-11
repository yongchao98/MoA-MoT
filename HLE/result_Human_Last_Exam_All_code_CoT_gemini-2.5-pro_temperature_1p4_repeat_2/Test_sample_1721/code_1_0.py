def solve_density_problem():
    """
    This script demonstrates the solution to find the largest density 'c'.
    The answer is c = 1/3, achieved by constructing the set A with numbers
    congruent to 1 modulo 3.
    """
    
    # We will demonstrate the construction for a sample size N.
    N = 1000

    # The optimal construction uses modulus m=3 and residue r=1.
    m = 3
    r = 1

    # Construct the set A = {n in {1,...,N} | n = r (mod m)}.
    # According to our plan, any sum of two elements from this set will be
    # congruent to 2 (mod 3), while squares can only be 0 or 1 (mod 3).
    # Thus, A+A will not contain any squares.
    A = [n for n in range(1, N + 1) if n % m == r]
    
    size_A = len(A)
    c_approx = size_A / N

    print(f"Demonstration for N = {N}:")
    print(f"The set A is constructed with numbers congruent to {r} modulo {m}.")
    print(f"The size of set A is |A| = {size_A}.")
    print(f"The calculated density for this N is |A|/N = {size_A}/{N} = {c_approx:.4f}")
    print("As N approaches infinity, this density converges to 1/3.")
    
    # The final answer is the theoretical limit, which is 1/3.
    # We print the components of this final fraction as requested.
    numerator = 1
    denominator = 3
    c_value = numerator / denominator

    print("\n-------------------------------------------")
    print("The final answer for the largest number c is based on the optimal construction.")
    print(f"The final equation for c is: {numerator} / {denominator}")
    print(f"This evaluates to: {c_value}")
    print("-------------------------------------------")

# Execute the function to show the result.
solve_density_problem()