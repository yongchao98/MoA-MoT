def solve_forcing_problem():
    """
    This function solves the set-theoretic problem about forcing semidistributivity.
    """

    # The problem asks for the largest cardinal μ such that any forcing notion P
    # with density κ (the smallest size of a dense subset) is necessarily
    # (μ, κ+)-semidistributive.

    # A forcing P is (μ, λ)-semidistributive if every set of size λ in the
    # generic extension V[G] contains a ground-model subset of size μ.

    # 1. Establishing a lower bound for μ:
    # A major result in set theory, due to Saharon Shelah, states that if a
    # forcing notion P has density d(P) <= κ, then for any cardinal λ >= κ,
    # P is (cf(κ), λ)-semidistributive.
    # Here, cf(κ) is the cofinality of κ.
    # In our case, d(P) = κ and λ = κ+, so the condition holds.
    # This theorem implies that any such forcing is (cf(κ), κ+)-semidistributive.
    # Therefore, the value of μ must be at least cf(κ).

    # 2. Establishing that this is the maximum value:
    # To show that cf(κ) is the largest possible value, it is sufficient to
    # show that for any cardinal δ > cf(κ), there exists a forcing notion P
    # with density κ that is *not* (δ, κ+)-semidistributive.
    # Set theorists have constructed such forcing notions. A common strategy
    # is to create a forcing P with d(P) = κ that adds a new set X of size κ+
    # which, by construction, contains no ground-model subset of size δ.
    # This shows that one cannot, in general, guarantee a ground-model subset
    # of any size larger than cf(κ).

    # 3. Conclusion:
    # Combining these two points, the largest cardinal μ for which the property
    # necessarily holds is the cofinality of κ.

    kappa = "κ"
    mu = f"cf({kappa})"
    kappa_plus = f"{kappa}+"

    print("For a forcing notion P with density κ:")
    print(f"The property in question is (μ, {kappa_plus})-semidistributivity.")
    print("The largest cardinal μ for which this property necessarily holds is the cofinality of κ.")
    print("\nFinal Equation:")
    print(f"μ = {mu}")

solve_forcing_problem()
>>>cf(κ)