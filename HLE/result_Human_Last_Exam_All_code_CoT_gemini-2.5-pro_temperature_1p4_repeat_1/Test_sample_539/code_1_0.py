import math

def check_options(k):
    """
    Checks the inequality ceil((k+1)/l) > k for the given options for l.
    """
    if not isinstance(k, int) or k <= 0:
        print(f"k must be a positive integer, but got {k}")
        return

    print(f"--- Checking for k = {k} ---")

    # Option A: l = k
    l_a = k
    val_a = math.ceil((k + 1) / l_a)
    holds_a = val_a > k
    print(f"Option A: l = {l_a}")
    print(f"Equation: ceil(({k} + 1) / {l_a}) > {k}")
    print(f"Result: ceil({(k + 1) / l_a:.2f}) > {k}  =>  {val_a} > {k}, which is {holds_a}")

    # Option B: l = k + 1
    l_b = k + 1
    val_b = math.ceil((k + 1) / l_b)
    holds_b = val_b > k
    print(f"\nOption B: l = {l_b}")
    print(f"Equation: ceil(({k} + 1) / {l_b}) > {k}")
    print(f"Result: ceil({(k + 1) / l_b:.2f}) > {k}  =>  {val_b} > {k}, which is {holds_b}")

    # Option C: l = k - 1
    if k - 1 > 0:
        l_c = k - 1
        val_c = math.ceil((k + 1) / l_c)
        holds_c = val_c > k
        print(f"\nOption C: l = {l_c}")
        print(f"Equation: ceil(({k} + 1) / {l_c}) > {k}")
        print(f"Result: ceil({(k + 1) / l_c:.2f}) > {k}  =>  {val_c} > {k}, which is {holds_c}")
    else:
        print("\nOption C: l = k - 1 is not a positive integer.")
    print("-" * 25)

# We test the inequality for k = 1, 2, 3
check_options(1)
check_options(2)
check_options(3)
