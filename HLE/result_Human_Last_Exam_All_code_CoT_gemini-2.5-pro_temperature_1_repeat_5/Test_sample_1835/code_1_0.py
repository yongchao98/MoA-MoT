def main():
    """
    This script models the logical relationship between understanding a singular
    proposition (Fa) and a universal one (∀x Fx).
    """

    # 1. Define the components of our logical world.
    # The predicate F is 'is_positive'.
    # The domain of discourse X is a set of numbers.
    # The individual 'a' is a number from that domain.
    predicate_F = lambda n: n > 0
    predicate_name = "is_positive"
    domain_X = {5, 10, 15, 20}
    individual_a = 5

    print(f"Let's assume the following:")
    print(f"Predicate F = {predicate_name}")
    print(f"Domain of individuals X = {domain_X}")
    print(f"A specific individual 'a' = {individual_a}\n")

    # 2. Model understanding the singular proposition 'Fa'.
    # If you understand 'Fa', you can apply the predicate F to the individual 'a'.
    result_fa = predicate_F(individual_a)
    print(f"--- Understanding 'Fa' ---")
    print(f"Evaluating the proposition F(a): {predicate_name}({individual_a})")
    print(f"Result: {result_fa}\n")

    # 3. Model the capacity to understand the universal proposition '∀x Fx'.
    # This requires having the concept of F (from 'Fa') and the universal quantifier '∀x'.
    # The script now combines these to check if F holds for ALL x in the domain X.
    print(f"--- Understanding '∀x Fx' (For all x, Fx) ---")
    print(f"To understand '∀x Fx', we must apply the predicate '{predicate_name}' to all individuals in the domain {domain_X}.")
    
    # Perform the universal check
    all_results = [predicate_F(x) for x in sorted(list(domain_X))]
    is_universally_true = all(all_results)
    
    print(f"The proposition '∀x Fx' is: {is_universally_true}\n")
    
    # 4. As requested, show the final logical "equation" with each number.
    # This expands the logic of the 'all()' function.
    print(f"--- The Final Equation for the Universal Check ---")
    boolean_checks = [f"{predicate_name}({x}) = {predicate_F(x)}" for x in sorted(list(domain_X))]
    full_equation = " AND ".join([f"{predicate_F(x)}" for x in sorted(list(domain_X))])
    
    print("The universal statement '∀x Fx' expands to:")
    # Print each individual check for clarity
    for check in boolean_checks:
        print(check)
    
    # Print the full equation showing each number (via its truth value)
    final_evaluation_string = " AND ".join(str(predicate_F(x)) for x in sorted(list(domain_X)))
    print(f"\nFinal Equation: {final_evaluation_string} = {is_universally_true}")


if __name__ == "__main__":
    main()