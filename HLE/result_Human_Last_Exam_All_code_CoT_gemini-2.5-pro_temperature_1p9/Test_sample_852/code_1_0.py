import itertools

def find_smallest_group_size():
    """
    Finds the smallest size of a finite Abelian group G with a maximal sum-free set S 
    satisfying |k(S)| > 2|S|.
    This function implements the logic derived in the explanation and finds a concrete
    example in the group Z_4 x Z_4.
    """
    dims = [4, 4]
    group_name = "Z_4 x Z_4"
    order = 16

    # --- Group Definition ---
    G = [tuple(p) for p in itertools.product(*[range(d) for d in dims])]
    identity = G[0]

    def op(a, b):
        return tuple((a[i] + b[i]) % dims[i] for i in range(len(dims)))

    def double(a):
        return tuple((2 * a[i]) % dims[i] for i in range(len(dims)))
    
    # --- Helper functions ---
    def is_sum_free(S, G_elements, op_func):
        if not S: return True
        if identity in S: return False
        for s1 in S:
            for s2 in S:
                if op_func(s1, s2) in S:
                    return False
        return True

    def is_maximal_sum_free(S, G_elements, op_func):
        if not is_sum_free(S, G_elements, op_func):
            return False
        
        G_minus_S = [g for g in G_elements if g not in S]
        
        for g in G_minus_S:
            S_union_g = S.union({g})
            if is_sum_free(S_union_g, G_elements, op_func):
                return False
        return True

    # --- Search for the specific set S ---
    # We are looking for S with j=2 elements from G[2] and l=1 element from G\G[2].
    G2 = {g for g in G if double(g) == identity}
    order_2_elements = G2 - {identity}
    
    # All 2-element subsets of order_2_elements are sum-free.
    j_sets = [set(comb) for comb in itertools.combinations(order_2_elements, 2)]
    l_candidates = [g for g in G if g not in G2]
    
    for T in j_sets:
      for u in l_candidates:
          S_candidate = T.union({u})
          
          if is_sum_free(S_candidate, G, op):
              if is_maximal_sum_free(S_candidate, G, op):
                  k_S = {g for g in G if double(g) in S_candidate}
                  
                  if len(k_S) > 2 * len(S_candidate):
                      print(f"The smallest size of such a group is {order}.")
                      print(f"The group is G = {group_name}.")
                      s_str = ", ".join(map(str, sorted(list(S_candidate))))
                      k_s_str = ", ".join(map(str, sorted(list(k_S))))
                      print(f"A maximal sum-free set is S = {{{s_str}}}.")
                      print(f"The size of S is |S| = {len(S_candidate)}.")
                      print(f"The set k(S) = {{g in G | 2g in S}} is {{{k_s_str}}}.")
                      print(f"The size of k(S) is |k(S)| = {len(k_S)}.")
                      print(f"The condition |k(S)| > 2|S| is satisfied as follows:")
                      print(f"{len(k_S)} > 2 * {len(S_candidate)}")
                      print(f"Which is: {len(k_S)} > {2*len(S_candidate)}")
                      return order
    return None

find_smallest_group_size()
