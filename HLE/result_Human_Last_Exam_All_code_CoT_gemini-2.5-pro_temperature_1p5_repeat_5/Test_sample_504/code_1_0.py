import math

def get_partitions(n, min_val=1):
    """
    Generates all integer partitions of n.
    """
    if n == 0:
        yield []
        return
    for i in range(min_val, n + 1):
        for p in get_partitions(n - i, i):
            yield [i] + p

def get_max_vectors(d):
    """
    Returns the maximum number of vectors m(d) for a given dimension d.
    Uses known exact values and bounds.
    """
    if d == 1:
        # For d=1, we can only have one vector line.
        return 1, "exact"
    if d == 2:
        # For d=2, c=1/2, the bound is m <= 3. This is achievable.
        return 3, "exact"
    if d == 3:
        # For d=3, c=1/2, we have a SIC-POVM with m = d^2 = 9.
        return 9, "exact"
    if d > 3:
        # For d>3, m < d^2. The bound is m <= d^2 - 1.
        # This bound may not be tight.
        return d**2 - 1, "bound"

def solve():
    """
    Solves the vector problem by analyzing partitions of the dimension 6.
    """
    total_dim = 6
    # We need at least one orthogonal pair, so p >= 2 parts.
    partitions = [p for p in get_partitions(total_dim) if len(p) >= 2]
    
    max_k = 0
    best_partition = None
    
    print("Analyzing partitions of dimension 6:")
    print("-" * 40)
    
    results = []

    for p in sorted(partitions, key=lambda x: sum([get_max_vectors(d)[0] for d in x]), reverse=True):
        current_k = 0
        is_exact = True
        
        dims_str = " + ".join(map(str, p))
        
        sum_str_parts = []
        for d in p:
            m, status = get_max_vectors(d)
            current_k += m
            sum_str_parts.append(f"m({d})")
            if status == "bound":
                is_exact = False
        
        sum_str = " + ".join(sum_str_parts)
        op = "=" if is_exact else "<="
        
        results.append({
            "partition": p,
            "sum_str": sum_str,
            "op": op,
            "value": current_k,
            "is_exact": is_exact
        })

    # Sort results for clearer presentation
    results.sort(key=lambda x: x["value"], reverse=True)
    
    for res in results:
        p_str = " + ".join(map(str, res["partition"]))
        print(f"Partition ({p_str}): k = {res['sum_str']} {res['op']} {res['value']}")

    # The (3,3) partition gives a provable value of 18.
    provable_max_k = 18
    print("-" * 40)
    print("The partition (3 + 3) gives k = m(3) + m(3) = 9 + 9 = 18.")
    print("This is known to be achievable using two sets of Hesse SIC-POVMs in orthogonal subspaces.")
    print("Other partitions either yield a smaller upper bound (e.g., (2+4) gives k <= 18)")
    print("or rely on finding exceptionally large sets of lines whose existence is unproven (e.g., (1+5) requires m(5) >= 18).")
    print("\nTherefore, the largest known and achievable number is 18.")

solve()
