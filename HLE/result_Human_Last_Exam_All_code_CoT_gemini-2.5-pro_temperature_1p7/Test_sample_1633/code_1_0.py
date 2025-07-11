import math

def solve():
    """
    Calculates the minimum fraction of rewiring operations m(n)/n needed
    to transform a network into an ultra-small world, and evaluates the options.
    """
    n = 1_000_000
    avg_k = 6.0
    min_k = 3
    # Per problem constraints, max degree is ceil(log(n))
    # log(1,000,000) = 13.8, so ceil is 14
    max_k = math.ceil(math.log(n))

    # A power-law distribution P(k) ~ k^-gamma is required for the
    # ultra-small world property L ~ loglog(n). We use gamma=2.5 as an example.
    gamma = 2.5
    
    # Define the degree range for the final graph
    degrees = list(range(min_k, max_k + 1))
    
    # To find the normalization constant C for P(k)=C*k^-gamma, we enforce that
    # the average degree must remain 6.
    # avg_k = sum(k * P(k)) = C * sum(k * k^-gamma)
    sum_k_pk_unnormalized = sum([k * (k**-gamma) for k in degrees])
    C = avg_k / sum_k_pk_unnormalized
    
    print("For n = {:,}".format(n))
    print(f"Degree range is from {min_k} to {max_k}.")
    print(f"Assuming a target power-law distribution P(k) ~ k^-{gamma}.")
    print(f"Calculated normalization constant C = {C:.4f}")

    # The final probability distribution for degrees
    p_k = {k: C * (k**-gamma) for k in degrees}
    
    # The transformation requires moving edge endpoints from nodes with degree < 6
    # to nodes with degree > 6. The total degree increase for high-degree nodes
    # is D_inc = n * sum_{k>6} (k - 6) * P(k).
    d_inc_per_node = sum([(k - avg_k) * p_k[k] for k in degrees if k > avg_k])
    
    # Each rewiring operation adds one edge, increasing the total degree sum by 2.
    # So, 2 * m(n) >= D_inc = n * d_inc_per_node.
    # This gives a lower bound on the ratio m(n)/n.
    m_div_n_lower_bound = d_inc_per_node / 2.0
    
    print("\n--- Calculating the minimum required rewiring operations ---")
    print(f"Total degree increase required for nodes with degree > 6 (per node): {d_inc_per_node:.4f}")
    
    equation_str = (f"m(n)/n >= (Sum_{{k=7}}^{{{max_k}}} (k - 6) * P(k)) / 2 "
                    f"= {d_inc_per_node:.4f} / 2 = {m_div_n_lower_bound:.4f}")
    print(f"The minimum number of rewirings m(n) must satisfy:")
    print(equation_str)
    
    # Statement H proposes m(n) >= n/6. Let's check our result.
    # n/6 is equivalent to m(n)/n >= 1/6
    h_bound = 1/6
    print(f"\nStatement H suggests m(n)/n >= 1/6, which is ~{h_bound:.4f}.")
    print(f"Our calculated lower bound is {m_div_n_lower_bound:.4f}.")
    if m_div_n_lower_bound >= h_bound:
        print("Our result is consistent with statement H.")
    else:
        print("Our result is not consistent with statement H.")

solve()
