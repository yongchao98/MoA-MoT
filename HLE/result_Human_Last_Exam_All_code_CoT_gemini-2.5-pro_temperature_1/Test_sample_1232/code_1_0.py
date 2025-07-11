def verify_potts_supermodularity_condition():
    """
    This function checks the supermodularity condition for the Potts model Hamiltonian,
    which is a sufficient condition for the positive correlations property.

    The condition is:
    I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d)) >= I(a=b) + I(c=d)
    where a,b,c,d are spin values from {1, ..., q}.

    We test a counterexample for q=3.
    """
    # Let q=3. The states are {1, 2, 3}.
    # Consider an edge (x,y) and two configurations, xi and eta.
    # Let xi(x) = a, xi(y) = b
    # Let eta(x) = c, eta(y) = d
    a = 2
    b = 2
    c = 1
    d = 3

    print("We check the supermodularity inequality for a counterexample with q=3.")
    print("The inequality required for the positive correlation property is:")
    print("I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d)) >= I(a=b) + I(c=d)\n")
    print(f"We choose the following values: a={a}, b={b}, c={c}, d={d}.\n")

    # Calculate the Right Hand Side (RHS) of the inequality
    # I(a=b) is 1 if a equals b, and 0 otherwise.
    # I(c=d) is 1 if c equals d, and 0 otherwise.
    I_ab = 1 if a == b else 0
    I_cd = 1 if c == d else 0
    rhs = I_ab + I_cd

    # Calculate the Left Hand Side (LHS) of the inequality
    max_ac = max(a, c)
    max_bd = max(b, d)
    I_max = 1 if max_ac == max_bd else 0

    min_ac = min(a, c)
    min_bd = min(b, d)
    I_min = 1 if min_ac == min_bd else 0
    lhs = I_max + I_min

    # Print the full equation with the calculated numbers
    print("Substituting these values into the inequality:")
    print(f"I(max({a},{c})=max({b},{d})) + I(min({a},{c})=min({b},{d})) >= I({a}={b}) + I({c}={d})")
    print(f"I({max_ac}={max_bd}) + I({min_ac}={min_bd}) >= {I_ab} + {I_cd}")
    print(f"       {I_max} + {I_min} >= {I_ab} + {I_cd}")
    final_inequality_str = f"       {lhs} >= {rhs}"
    print(final_inequality_str)

    if lhs >= rhs:
        print("\nThe inequality holds for this case.")
    else:
        print("\nThe inequality FAILS for this case.")
        print("\nThis failure for q=3 implies the positive correlation property does not hold in general for graphs with at least one edge (i.e., max degree >= 1).")
        print("Since the property must hold for ALL q >= 2, the largest possible max degree 'd' must be 0.")

verify_potts_supermodularity_condition()