def explain_mf_dimensions():
    """
    Explains the dimensional structure of Type-1, Type-2, and Type-3
    fuzzy membership functions by analyzing their domains.
    """

    # The domain of a function defines its dimensional input structure.
    # Let's represent the variables in the domain for each MF type.
    # x: The primary variable from the universe of discourse.
    # u: The secondary variable, for the membership of x.
    # v: The tertiary variable, for the membership of u.

    domain_variables_type1 = ["x"]
    domain_variables_type2 = ["x", "u"]
    domain_variables_type3 = ["x", "u", "v"]

    num_vars_type1 = len(domain_variables_type1)
    num_vars_type2 = len(domain_variables_type2)
    num_vars_type3 = len(domain_variables_type3)

    print("--- Analysis of Fuzzy Membership Function Domains ---")
    print(f"Type-1 MF domain is a function of {num_vars_type1} variable: f({', '.join(domain_variables_type1)})")
    print(f"Type-2 MF domain is a function of {num_vars_type2} variables: f({', '.join(domain_variables_type2)})")
    print(f"Type-3 MF domain is a function of {num_vars_type3} variables: f({', '.join(domain_variables_type3)})")
    print("-" * 51)

    print("\nThe fundamental difference in dimensional structure from Type-2 to Type-3 is the change in the function's domain.")
    print(f"It expands from a {num_vars_type2}-variable domain to a {num_vars_type3}-variable domain.")
    print("\nThis directly corresponds to the following answer choice:")
    print("C. Expanded to three-variable domain")

if __name__ == "__main__":
    explain_mf_dimensions()