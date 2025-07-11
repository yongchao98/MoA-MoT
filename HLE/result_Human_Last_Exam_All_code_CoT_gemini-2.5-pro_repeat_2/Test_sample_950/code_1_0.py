import sys

def solve_grassmannian_cohomology_rank():
    """
    This script calculates the rank of the torsion subgroup of the integral
    cohomology ring of the space of 3-subspaces of R^5.
    """

    # 1. Define the space
    k = 3
    n = 5
    space_name = f"Gr({k}, {n})"
    space_desc = f"the space of {k}-subspaces of R^{n}"

    print(f"The problem asks for the rank of the torsion subgroup of the integral cohomology ring of {space_desc}.")
    print(f"This space is the real Grassmannian manifold {space_name}.")
    print("-" * 30)

    # 2. Explain the theory
    print("Step 1: Understand the structure of the integral cohomology ring H*(Gr(k, n); Z).")
    print("From algebraic topology, it is known that the integral cohomology of any real Grassmannian Gr(k, n) is torsion-free.")
    print("This is because the Grassmannian has a CW decomposition (the Schubert cell decomposition) for which the cellular chain complex has trivial boundary maps.")
    print("This implies that the integral homology groups H_i(Gr(k, n); Z) are free abelian groups.")
    print("By the Universal Coefficient Theorem, the integral cohomology groups H^i(Gr(k, n); Z) are also free abelian groups.")
    print("\nThis means that for each degree i, H^i(Gr(k, n); Z) is isomorphic to a direct sum of copies of Z, i.e., Z^b_i, where b_i is the i-th Betti number.")
    print("-" * 30)

    # 3. Calculate the Betti numbers for Gr(3, 5) as a demonstration
    print(f"Step 2: Let's calculate the Betti numbers for {space_name} to see the structure explicitly.")
    print("The Betti numbers are determined by counting Schubert symbols, which are partitions fitting into a k x (n-k) grid.")
    print(f"For {space_name}, k={k} and n={n}, so we are looking for partitions into at most {k} parts, with each part at most {n-k}.")

    def generate_schubert_symbols(num_parts, max_val):
        """
        Generates sequences (s_1, ..., s_k) with max_val >= s_1 >= ... >= s_k >= 0.
        """
        if num_parts == 0:
            yield []
            return
        
        # To avoid generating duplicates, we ensure the next part is <= the current part.
        # The parameter `current_max` enforces this.
        def helper(p, current_max):
            if p == 0:
                yield []
                return
            for val in range(current_max, -1, -1):
                for sub_symbol in helper(p - 1, val):
                    yield [val] + sub_symbol
        
        yield from helper(num_parts, max_val)

    dim = k * (n - k)
    betti_numbers = [0] * (dim + 1)
    max_part_val = n - k
    
    for symbol in generate_schubert_symbols(k, max_part_val):
        symbol_dim = sum(symbol)
        if symbol_dim <= dim:
            betti_numbers[symbol_dim] += 1
            
    for i, b_i in enumerate(betti_numbers):
        if b_i > 0:
            print(f"H^{i}({space_name}; Z) = Z^{b_i}")

    print("-" * 30)

    # 4. Determine the torsion subgroup and its rank
    print("Step 3: Analyze the torsion subgroup and its rank.")
    print(f"The entire cohomology ring H*({space_name}; Z) is the direct sum of these free abelian groups.")
    print("A direct sum of torsion-free groups is itself torsion-free.")
    print("Therefore, the torsion subgroup of H*({space_name}; Z) is the trivial group {0}.")
    print("\nThe rank of an abelian group is the number of independent generators of its free part. For any torsion group (a group where every element has finite order), the rank is 0.")
    print("Since the torsion subgroup is {0}, its rank is 0.")
    print("-" * 30)

    # 5. State the final answer
    rank_of_torsion_subgroup = 0
    print("Final Answer:")
    print(f"The rank of the torsion subgroup of the integral cohomology ring of the space of {k}-subspaces of R^{n} is {rank_of_torsion_subgroup}.")

if __name__ == '__main__':
    solve_grassmannian_cohomology_rank()