import math

def solve_for_K():
    """
    This function solves for the value of K based on the asymptotic behavior of the slice rank.
    
    The exact slice rank is known to be 3^n - 2^n + 1.
    For large n, this is dominated by 3^n.
    The given formula for the slice rank is (3 / 2^K)^n * e^(o(n)).
    For large n, this is dominated by (3 / 2^K)^n.
    
    Equating the dominant terms, we get:
    3^n = (3 / 2^K)^n
    
    Taking the n-th root of both sides:
    3 = 3 / 2^K
    
    Solving for K:
    2^K = 3 / 3
    2^K = 1
    K = log2(1)
    """
    
    # Step 1: State the asymptotic equivalence
    print("The slice rank of the tensor is known to be 3^n - 2^n + 1.")
    print("For large n, the asymptotic behavior is dominated by the largest term, 3^n.")
    print("Slice Rank ~ 3^n")
    print("")
    print("The given slice rank formula is (3/2^K)^n * e^(o(n)).")
    print("The asymptotic behavior of this formula is dominated by (3/2^K)^n.")
    print("Slice Rank ~ (3/2^K)^n")
    print("")
    
    # Step 2: Equate the asymptotic expressions
    print("Equating the asymptotic forms:")
    print("3^n = (3 / 2^K)^n")
    print("")

    # Step 3: Simplify the equation
    print("Taking the n-th root of both sides:")
    base_exact = 3
    print(f"{base_exact} = 3 / (2^K)")
    print("")

    # Step 4: Solve for 2^K
    print("Solving for 2^K:")
    val_2k = base_exact / base_exact
    print(f"2^K = {base_exact} / {base_exact}")
    print(f"2^K = {val_2k}")
    print("")

    # Step 5: Solve for K
    print("Solving for K by taking log base 2:")
    # Using math.log2, though for 1 it's trivial.
    K = math.log2(val_2k)
    print(f"K = log2({val_2k})")
    print(f"K = {K}")

solve_for_K()