def main():
    """
    This script explains why Functional Extensionality is inconsistent
    with the given recursion scheme in dependent type theory.
    """

    # Let's use strings to represent terms in our logical system.
    # Consider an arbitrary function 'f'.
    f = "f"

    # 1. The Flawed Recursion Rule
    # The subterm rule allows defining a recursive function `Omega` such that:
    # Omega(f) is defined as a recursive call on the eta-expansion of f.
    # Let's represent the eta-expansion of `f` as `(lambda x. f(x))`.
    f_eta = f"(lambda x. {f}(x))"

    print("The system allows a non-terminating function `Omega` to be defined.")
    print("Its computational rule is a one-way reduction:")
    reduction_rule = f"Omega({f}) --> Omega({f_eta})"
    print(reduction_rule)
    print("-" * 40)

    # 2. The Axiom of Functional Extensionality
    # This axiom states that if two functions behave identically, they are equal.
    # A function `f` and its eta-expansion `f_eta` behave identically by definition.
    # Therefore, functional extensionality allows us to prove their equality.
    print("Functional Extensionality lets us prove the following theorem:")
    equality_from_funext = f"{f} = {f_eta}"
    print(equality_from_funext)
    print("-" * 40)

    # 3. The Contradiction
    # In a consistent logic, if we have a proof that `t = u`, then `t` and `u`
    # must be interchangeable. We can apply `Omega` to both sides of our proven equality.
    # This leads to a new provable equality:
    paradox_equation_lhs = f"Omega({f})"
    paradox_equation_rhs = f"Omega({f_eta})"
    
    print("Combining the rule and the axiom leads to a paradoxical provable equality:")
    print(f"{paradox_equation_lhs} = {paradox_equation_rhs}")
    print("\nThis is a contradiction because the term on the left reduces to the term on the right.")
    print("A system that allows proving `t = u` when `t` reduces to `u` is inconsistent,")
    print("as it can be used to construct a proof for any proposition, including False.")


if __name__ == "__main__":
    main()
