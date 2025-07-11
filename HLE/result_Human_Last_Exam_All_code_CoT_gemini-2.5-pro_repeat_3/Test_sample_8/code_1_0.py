import sys

def solve_conormal_space():
    """
    This script determines the conormal space of R(sigma)f based on the principle of elliptic regularity.
    """
    # 1. Define the order of the differential operator
    # The operator is P = Box_g - sigma(sigma-1). Box_g is the wave operator, which is a
    # second-order differential operator. Thus, the order of P is 2.
    operator_order = 2

    # 2. Define the regularity of the source term f
    # f is in A^{2+alpha}(X), so its conormal regularity order is 2+alpha.
    # We represent this symbolically.
    source_regularity_val = 2
    source_regularity_sym = "alpha"

    # 3. Apply the principle of elliptic regularity
    # The regularity of the solution u = R(sigma)f is the sum of the
    # regularity of the source f and the order of the elliptic operator.
    # New regularity = (source regularity) + (operator order)
    solution_regularity_val = source_regularity_val + operator_order

    # 4. Print the derivation
    print("Derivation of the conormal space for R(sigma)f:")
    print("-" * 50)
    print("The problem is to find the conormal space of u = R(sigma)f, which is the solution to the equation:")
    print("  (Box_g - sigma(sigma-1))u = f\n")

    print(f"1. The differential operator P = Box_g - sigma(sigma-1) is elliptic of order m.")
    print(f"   m = {operator_order}\n")

    print(f"2. The source term f belongs to the conormal space A^(s)(X), where s is the order of regularity.")
    print(f"   s = {source_regularity_val} + {source_regularity_sym}\n")

    print("3. According to the principle of elliptic regularity, the solution u = R(sigma)f has a regularity order s' = s + m.")
    print("   This means the resolvent R(sigma), being the inverse of P, is an operator of order -m = -2,")
    print("   and it increases the regularity of the function it is applied to by 2.\n")

    print("4. The final regularity order s' is calculated as follows:")
    print(f"   s' = s + m")
    print(f"   s' = ({source_regularity_val} + {source_regularity_sym}) + {operator_order}")
    print(f"   s' = {solution_regularity_val} + {source_regularity_sym}\n")

    print("Conclusion:")
    print(f"Therefore, R(sigma)f belongs to the conormal space A^({solution_regularity_val}+{source_regularity_sym})(X).")

if __name__ == "__main__":
    solve_conormal_space()
