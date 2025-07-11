import itertools

def solve_module_problem():
    """
    Finds the unique tau-tilting module that is not a slice for the
    path algebra A = C(1->2->3).
    """

    # The 6 indecomposable modules over A
    indecs = ["S1", "S2", "S3", "M12", "M23", "M13"]

    # The support of each indecomposable module, indicating which simple
    # modules are involved in its structure.
    supports = {
        "S1": {1},
        "S2": {2},
        "S3": {3},
        "M12": {1, 2},
        "M23": {2, 3},
        "M13": {1, 2, 3}
    }

    # The Auslander-Reiten translation (tau) acts on the indecomposables.
    # For this algebra, the only non-zero translation is tau(S2) = S1.
    # tau(M) = 0 for all other indecomposable modules M since they are
    # either projective or injective.
    tau_map = {"S2": "S1"}

    # To check for tau-rigidity, Hom(T, tau(T)) = 0, we need to know
    # which modules X have Hom(X, S1) != 0.
    # These are the modules that have S1 as a quotient.
    hom_to_S1_non_zero = {"S1", "M12", "M13"}

    # Dimension vectors for each indecomposable module.
    dim_vectors = {
        "S1": (1, 0, 0), "S2": (0, 1, 0), "S3": (0, 0, 1),
        "M12": (1, 1, 0), "M23": (0, 1, 1), "M13": (1, 1, 1)
    }

    def is_tau_rigid(module_summands):
        """A module is tau-rigid if Hom(M, tau(M)) = 0."""
        # This condition fails only if S2 is a summand, and another summand X
        # has Hom(X, tau(S2)) = Hom(X, S1) != 0.
        if "S2" in module_summands:
            for m in module_summands:
                if m in hom_to_S1_non_zero:
                    return False
        return True

    def is_sincere(module_summands):
        """A module is sincere if its support covers all simple modules {1, 2, 3}."""
        total_support = set()
        for m in module_summands:
            total_support.update(supports[m])
        return total_support == {1, 2, 3}

    # The number of vertices in the quiver.
    n = 3
    
    non_slice_modules = []

    # Iterate through all combinations of n indecomposable modules.
    for combo in itertools.combinations(indecs, n):
        module_summands = set(combo)
        # A module is tau-tilting if it is tau-rigid and has n summands.
        if is_tau_rigid(module_summands):
            # A slice module must be sincere. We are looking for one that is not.
            if not is_sincere(module_summands):
                non_slice_modules.append(sorted(list(module_summands)))

    # The problem implies there is a unique such module.
    if len(non_slice_modules) == 1:
        result_module = non_slice_modules[0]
        
        print("The unique tau-tilting module that is not a slice is T, where:")
        # The module is a direct sum of its indecomposable summands.
        module_string = " \u2295 ".join(result_module)
        print(f"T = {module_string}")

        # Map to other common names: S3 is the projective P3, M23 is P2.
        pretty_names = {"S2": "S_2", "S3": "P_3", "M23": "P_2"}
        pretty_string = " \u2295 ".join([pretty_names.get(m, m) for m in result_module])
        print(f"Also known as: T = {pretty_string}")
        
        print("\nThe equation for its dimension vector is:")
        
        total_dim_vec = [0, 0, 0]
        equation_parts = []
        for m in result_module:
            vec = dim_vectors[m]
            total_dim_vec[0] += vec[0]
            total_dim_vec[1] += vec[1]
            total_dim_vec[2] += vec[2]
            equation_parts.append(str(vec))

        # Outputting each number in the final equation, as requested.
        sum_equation = " + ".join(equation_parts)
        result_vec_str = f"({total_dim_vec[0]}, {total_dim_vec[1]}, {total_dim_vec[2]})"
        print(f"{sum_equation} = {result_vec_str}")
        
    else:
        print(f"Error: Found {len(non_slice_modules)} modules, but expected 1.")
        if non_slice_modules:
            print("Found modules:", non_slice_modules)

solve_module_problem()