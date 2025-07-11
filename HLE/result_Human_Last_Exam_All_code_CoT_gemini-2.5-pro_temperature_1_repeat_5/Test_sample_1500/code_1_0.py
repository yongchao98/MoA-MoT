def analyze_betti_number_for_sun(n, multiplicities):
    """
    Analyzes the second Betti number for a coadjoint orbit of SU(n).

    The coadjoint orbit is determined by an element lambda, whose eigenvalue
    multiplicities are given. The second Betti number b2 is k-1, where k
    is the number of distinct eigenvalues. The question claims b2 is always n-1.

    Args:
        n (int): The dimension for SU(n).
        multiplicities (list of int): A list representing the multiplicities
                                      of the eigenvalues of lambda.
    """
    if sum(multiplicities) != n:
        raise ValueError("Sum of multiplicities must be equal to n.")

    # The number of distinct eigenvalues, k.
    k = len(multiplicities)

    # For a coadjoint orbit O_lambda = SU(n)/P, where lambda has k
    # distinct eigenvalues, the second Betti number b2 is k-1.
    b2 = k - 1

    # The rank of SU(n) is n-1. This is the value claimed by the question.
    claimed_b2 = n - 1

    print(f"--- Analysis for SU({n}) with eigenvalue multiplicities {multiplicities} ---")
    print(f"The element lambda has k = {k} distinct eigenvalues.")
    print(f"The true second Betti number is b2 = k - 1 = {k} - 1 = {b2}.")
    print(f"The number claimed in the question is n - 1 = {n} - 1 = {claimed_b2}.")

    if b2 == claimed_b2:
        print(f"In this case ({b2} == {claimed_b2}), the statement holds (lambda is regular).")
    else:
        print(f"Since {b2} != {claimed_b2}, this case provides a counterexample to the statement.")
        print("The statement 'b2 is always n-1' is therefore false.")

if __name__ == '__main__':
    print("Analyzing the statement: 'For G = SU(n), b2(O_lambda) is always n-1'.")
    print("We test this with a counterexample where lambda is singular.")
    
    # Counterexample: G = SU(3), lambda is singular.
    # This corresponds to an orbit like the Grassmannian Gr(2, 3), which is CP^2.
    # The multiplicities of eigenvalues are (2, 1).
    n_counterexample = 3
    multiplicities_counterexample = [2, 1]
    
    analyze_betti_number_for_sun(n_counterexample, multiplicities_counterexample)
