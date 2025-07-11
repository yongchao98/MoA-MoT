import math

class Infinity:
    """A simple class to represent infinity for arithmetic operations."""
    def __add__(self, other):
        return self

    def __radd__(self, other):
        return self

    def __mul__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return 0
        return self

    def __rmul__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return 0
        return self

    def __str__(self):
        return "inf"

    def __repr__(self):
        return "inf"

inf = Infinity()

def get_dim_Hb_T(k):
    """
    Returns the dimension of the k-th bounded cohomology group H_b^k(T; R)
    of Thompson's group T with trivial real coefficients.
    """
    if not isinstance(k, int) or k < 0:
        raise ValueError("Degree k must be a non-negative integer.")
    if k == 0:
        return 1
    if k == 1:
        return 0
    # For k >= 2, the dimension is infinite.
    return inf

def compute_dimension():
    """
    Computes the dimension of H_b^4(T x T; R) using the Künneth formula.
    """
    n = 4
    
    print("To compute the dimension of the degree 4 bounded cohomology group of T x T,")
    print("we use the Künneth formula for bounded cohomology.")
    print("T is a non-amenable group, so the formula applies.")
    print()
    print(f"dim H_b^{n}(T x T) = sum over p+q={n} of [ dim(H_b^p(T)) * dim(H_b^q(T)) ]")
    print("-" * 50)
    
    total_dim = 0
    
    # Store components for printing the final equations
    symbolic_terms = []
    numeric_terms = []
    evaluated_terms = []

    for p in range(n + 1):
        q = n - p
        dim_p = get_dim_Hb_T(p)
        dim_q = get_dim_Hb_T(q)
        
        term_dim = dim_p * dim_q
        
        symbolic_terms.append(f"dim(H_b^{p}(T))*dim(H_b^{q}(T))")
        numeric_terms.append(f"{dim_p}*{dim_q}")
        evaluated_terms.append(f"{term_dim}")
        
        total_dim += term_dim
        
    print("The sum expands to:")
    print(" + ".join(symbolic_terms))
    print("\nSubstituting the known dimensions:")
    print("= " + " + ".join(numeric_terms))
    print("\nEvaluating each term:")
    print("= " + " + ".join(evaluated_terms))
    print("\nSumming the terms:")
    print(f"= {total_dim}")
    print("-" * 50)
    print(f"The dimension of the degree 4 bounded cohomology group of T x T is {total_dim}.")

compute_dimension()