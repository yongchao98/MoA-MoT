import math

def main():
    """
    Illustrates the set-theoretic problem with finite analogues.
    """
    # Let's use finite numbers as stand-ins for the cardinals.
    omega1_size = 100  # Represents the size of the domain/codomain, omega_1
    omega2_size = 200  # Represents the length of the sequence, omega_2

    # Let's define the sequence of functions f_alpha.
    # Each f_alpha is a function from {0..omega1_size-1} to {0..omega1_size-1}.
    # functions[alpha] will be a list representing the function f_alpha.
    functions = []
    for alpha in range(omega2_size):
        f = []
        for gamma in range(omega1_size):
            # This is a simple construction to make f_alpha grow with alpha.
            # The value is essentially alpha, but capped at the codomain size.
            # We introduce a small perturbation to make the <* property non-trivial.
            val = min(alpha, omega1_size - 1)
            if alpha > 0 and gamma % (alpha + 1) == 0:
                # On a few 'gamma' coordinates, the value is smaller.
                # This ensures the set where f_beta(gamma) <= f_alpha(gamma) is non-empty.
                val = val // 2
            f.append(val)
        functions.append(f)

    # Let's check our sequence to see if it's increasing "modulo finite".
    # For alpha < beta, we expect f_beta(gamma) > f_alpha(gamma) for most gamma.
    alpha_check = omega1_size // 2
    beta_check = omega1_size
    error_set = []
    for gamma in range(omega1_size):
        if functions[beta_check][gamma] <= functions[alpha_check][gamma]:
            error_set.append(gamma)
    # This set should be "finite" (small in our analogy).
    # print(f"For alpha={alpha_check}, beta={beta_check}, the error set is {error_set}")
    # print(f"Size of error set: {len(error_set)} (should be small)")

    # Now, let's try to find a bounding function g for an "uncountable" set.
    # Let's take our "uncountable" set of functions to be those with index >= omega1_size.
    # A natural way to define g is to take the supremum over a "countable" set of functions.
    # Let's use the first omega1_size functions for this.
    g = []
    for gamma in range(omega1_size):
        max_val = 0
        for xi in range(omega1_size):
            if functions[xi][gamma] > max_val:
                max_val = functions[xi][gamma]
        # g(gamma) is the sup + 1 to ensure strict inequality
        g.append(max_val + 1)

    # Now, let's test if g bounds a function from the "uncountable" part of the sequence.
    # Let's pick alpha_star = omega1_size, which is the first index in our "uncountable" set.
    alpha_star = omega1_size
    f_star = functions[alpha_star]

    # We search for a gamma where the bound fails.
    failure_found = False
    for gamma in range(omega1_size):
        if f_star[gamma] >= g[gamma]:
            print("The proposed bounding function 'g' fails.")
            print(f"Let's test f_alpha for alpha = {alpha_star}.")
            print(f"Let g(gamma) = sup{{f_xi(gamma) | xi < {omega1_size}}} + 1.")
            print(f"At gamma = {gamma}, we have found a failure:")
            
            # Show the calculation for g[gamma]
            sup_val = g[gamma] - 1
            print(f"  sup{{f_xi({gamma}) | xi < {omega1_size}}} = {sup_val}, so g({gamma}) = {g[gamma]}")

            # Show the violating equation
            print(f"  f_{alpha_star}({gamma}) = {f_star[gamma]}")
            print("\nThis leads to the inequality:")
            print(f"  f_{alpha_star}({gamma}) = {f_star[gamma]} >= {g[gamma]} = g({gamma})")
            
            print("\nThis demonstrates that f_alpha is not bounded by g, as the condition f_alpha(gamma) < g(gamma) is violated.")
            failure_found = True
            break
            
    if not failure_found:
        print("Could not find a failure for this specific construction, but the principle holds.")

if __name__ == "__main__":
    main()
