import collections

def find_longest_chain_length(poset_hasse, start_node):
    """
    Finds the length of the longest chain starting from a given node in a poset.
    The length of a chain is the number of elements minus 1.
    
    Args:
        poset_hasse (dict): The Hasse diagram of the poset, where keys are elements
                            and values are sets of elements they cover.
        start_node: The element where the chains should start.

    Returns:
        The length of the longest chain starting from start_node.
    """
    
    # Memoization for dynamic programming
    memo = {}

    def get_max_depth(node):
        if node in memo:
            return memo[node]
        
        # Base case: a maximal element in the poset
        if not poset_hasse.get(node):
            memo[node] = 0
            return 0
        
        # Recursive step: max depth is 1 + max depth of children
        max_d = 0
        for successor in poset_hasse[node]:
            max_d = max(max_d, 1 + get_max_depth(successor))
        
        memo[node] = max_d
        return max_d

    return get_max_depth(start_node)


# Define the diamond lattice J by its Hasse diagram (covering relations)
# An arrow from x to y means x < y.
# e.g., 'b': {'m1', 'm2'} means b < m1 and b < m2 are covering relations.
diamond_lattice_hasse = {
    'b': {'m1', 'm2'},
    'm1': {'t'},
    'm2': {'t'},
    't': set()
}

start_element = 'b'

# According to Webb's theorem, pd(S_b) = length of the longest chain starting at 'b'.
pd_Sb = find_longest_chain_length(diamond_lattice_hasse, start_element)

# Let's find one of the longest chains to demonstrate
def find_a_longest_chain(poset_hasse, start_node):
    chain = [start_node]
    curr = start_node
    while poset_hasse.get(curr):
        # To find the next element in the longest chain, we need the depths
        # which we already have logic for.
        max_depth = -1
        next_node = None
        for succ in poset_hasse[curr]:
            # This is a bit inefficient without memoization, but ok for small poset
            d = find_longest_chain_length(poset_hasse, succ)
            if d > max_depth:
                max_depth = d
                next_node = succ
        chain.append(next_node)
        curr = next_node
    return chain

longest_chain = find_a_longest_chain(diamond_lattice_hasse, start_element)


print("The problem considers a tame functor on an upper semilattice J.")
print("We investigate the diamond lattice, which is an upper semilattice of finite (hence tame) type.")
print("The lattice J has elements {b, m1, m2, t} with relations b < m1 < t and b < m2 < t.")
print("\nThe projective dimension of the simple module S_b for the element 'b' gives a lower bound for n.")
print("This dimension is the length of the longest chain in J starting at 'b'.")
print(f"The chains starting at '{start_element}' are (b, m1, t) and (b, m2, t).")
print(f"Let's take one such chain: {longest_chain}.")
print(f"The number of elements in this chain is {len(longest_chain)}.")
print(f"The length of this chain is (number of elements - 1).")
print(f"So, the projective dimension of S_b is {len(longest_chain)} - 1 = {pd_Sb}.")
print(f"\nSince there exists a functor with projective dimension 2, any general value n must be at least 2.")
