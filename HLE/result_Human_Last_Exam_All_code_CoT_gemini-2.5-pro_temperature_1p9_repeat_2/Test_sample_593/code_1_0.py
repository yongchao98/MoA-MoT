import sys

def solve():
    """
    This function takes three integer inputs t_H, t_G, and k,
    and calculates the tight upper bound on the treewidth of the combined graph F.
    """
    try:
        # Read inputs for t_H, t_G, and k from stdin
        # For demonstration, we'll use example values.
        # To run with user input:
        # inputs = sys.stdin.readline().split()
        # t_H = int(inputs[0])
        # t_G = int(inputs[1])
        # k = int(inputs[2])

        # Example values from the reasoning
        # Case 1: K_k,k example, where k=5
        # t_H = 1, t_G = 4, k = 5. Correct answer should be 5.
        t_H = 1
        t_G = 4
        k = 5
        
        # Calculate the bound for the first example
        # Bound is max(t_H, t_G, k-1) + 1
        # The term k-1 is meaningful only if k > 0. If k=0, there are no shared vertices,
        # and the treewidth is simply max(t_H, t_G). The formula still works if we
        # consider max with a negative number.
        if k == 0:
            bound_val_1 = max(t_H, t_G)
        else:
            bound_val_1 = max(t_H, t_G, k - 1) + 1
            
        print("--- Example 1: F = K_k,k ---")
        print(f"Let treewidth of H be t_H = {t_H}")
        print(f"Let treewidth of G be t_G = {t_G}")
        print(f"Let the number of shared vertices be k = {k}")
        print("The treewidth of the combined graph F is upper-bounded by max(t_H, t_G, k-1) + 1.")
        print(f"Upper bound = max({t_H}, {t_G}, {k}-1) + 1 = max({t_H}, {t_G}, {k-1}) + 1 = {max(t_H, t_G, k-1)} + 1 = {bound_val_1}")
        print("-" * 20)

        # Case 2: K_10 and P_3 example
        # t_H = 9, t_G = 1, k = 2. Correct answer should be 9. The bound is 10.
        t_H = 9
        t_G = 1
        k = 2

        # Calculate the bound for the second example
        if k == 0:
            bound_val_2 = max(t_H, t_G)
        else:
            bound_val_2 = max(t_H, t_G, k - 1) + 1
            
        print("--- Example 2: F is union of K_10 and P_3 ---")
        print(f"Let treewidth of H be t_H = {t_H}")
        print(f"Let treewidth of G be t_G = {t_G}")
        print(f"Let the number of shared vertices be k = {k}")
        print("The treewidth of the combined graph F is upper-bounded by max(t_H, t_G, k-1) + 1.")
        print(f"Upper bound = max({t_H}, {t_G}, {k}-1) + 1 = max({t_H}, {t_G}, {k-1}) + 1 = {max(t_H, t_G, k-1)} + 1 = {bound_val_2}")
        print("-" * 20)


    except (ValueError, IndexError):
        print("Please provide three integers for t_H, t_G, and k.", file=sys.stderr)

solve()