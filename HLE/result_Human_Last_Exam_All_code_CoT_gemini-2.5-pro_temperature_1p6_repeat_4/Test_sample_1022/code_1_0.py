import sys

def main():
    """
    Computes the dimension of the degree 4 bounded cohomology group of T x T.
    """
    degree = 4

    print("This script computes the dimension of the degree 4 bounded cohomology group")
    print("of the group G = T x T, where T is Thompson's group T, with trivial real coefficients.")
    print("The quantity to compute is dim H_b^4(T x T; R).\n")

    print("Step 1: The Künneth Formula for Bounded Cohomology")
    print("-------------------------------------------------------")
    print("For two groups G1 and G2, the Künneth formula provides an isomorphism:")
    print("H_b^n(G1 x G2; R) ~= SUM_{p+q=n} (H_b^p(G1; R) (x) H_b^q(G2; R))")
    print("Taking dimensions, we get:")
    print("dim H_b^n(G1 x G2; R) = SUM_{p+q=n} dim H_b^p(G1; R) * dim H_b^q(G2; R)\n")
    
    print("Step 2: Bounded Cohomology of Thompson's Group T")
    print("-------------------------------------------------")
    print("The dimensions of the bounded cohomology groups of T with real coefficients are known:")
    
    # This function provides the known dimensions of H_b^n(T; R)
    def get_dim_Hb_T(n):
        if n == 0:
            return 1
        if n == 2:
            return 1
        # For n=1 and n>=3, the dimension is 0
        return 0

    print(f"dim H_b^0(T; R) = {get_dim_Hb_T(0)}")
    print(f"dim H_b^1(T; R) = {get_dim_Hb_T(1)}")
    print(f"dim H_b^2(T; R) = {get_dim_Hb_T(2)}")
    print("dim H_b^n(T; R) = 0 for n >= 3\n")

    print(f"Step 3: Applying the Künneth Formula for n = {degree}")
    print("-------------------------------------------------")
    print("We want to compute dim H_b^4(T x T; R). Using the formula with G1=T, G2=T, n=4:")
    print("dim H_b^4(T x T) = "
          "dim H_b^0(T)*dim H_b^4(T) + "
          "dim H_b^1(T)*dim H_b^3(T) + "
          "dim H_b^2(T)*dim H_b^2(T) + "
          "dim H_b^3(T)*dim H_b^1(T) + "
          "dim H_b^4(T)*dim H_b^0(T)\n")

    print("Step 4: Calculation")
    print("-------------------")
    
    terms_values = []
    # Loop through all p from 0 to degree
    for p in range(degree + 1):
        q = degree - p
        dim_p = get_dim_Hb_T(p)
        dim_q = get_dim_Hb_T(q)
        term_value = dim_p * dim_q
        terms_values.append(term_value)
        
        print(f"Term for (p={p}, q={q}): dim H_b^{p}(T) * dim H_b^{q}(T) = {dim_p} * {dim_q} = {term_value}")

    total_dimension = sum(terms_values)
    
    # Construct the final equation string with all numbers
    final_equation_str = " + ".join(map(str, terms_values))
    
    print("\nSumming the terms:")
    print(f"dim H_b^4(T x T; R) = {final_equation_str} = {total_dimension}")
    
    print("\nFinal Result:")
    print("--------------")
    print(f"The dimension of the degree 4 bounded cohomology group of T x T is {total_dimension}.")
    
    # Output the final answer in the specified format
    # Redirecting to stderr to not interfere with final answer format.
    print("\n<<<1>>>", file=sys.stderr)


if __name__ == "__main__":
    main()