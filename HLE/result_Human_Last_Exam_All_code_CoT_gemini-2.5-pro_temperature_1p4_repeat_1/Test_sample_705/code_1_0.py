import itertools

def demonstrate_exponential_blowup():
    """
    This function demonstrates that the relaxed simulation can have exponentially
    growing computational cost compared to an ordinary simulation. We use the
    example of an n-bit counter.
    """
    n = 10  # Number of components (bits in our counter)

    # Let each V_k be the set of possible values for the k-th bit.
    # To make them disjoint, we use tuples: V_k = {(k, 0), (k, 1)}.
    V = {k: {(k, 0), (k, 1)} for k in range(n)}
    
    # The full domain D is the union of all V_k
    D_domain = set().union(*V.values())

    def to_int(s):
        """Converts a state tuple to an integer."""
        num = 0
        for i in range(n):
            bit_val = s[i][1]
            num += bit_val * (2**i)
        return num

    def from_int(num):
        """Converts an integer to a state tuple."""
        s = []
        for i in range(n):
            bit_val = 1 if (num >> i) & 1 else 0
            s.append((i, bit_val))
        return tuple(s)

    def f(s):
        """The simulator function: increments the integer value of the state."""
        current_val = to_int(s)
        next_val = (current_val + 1) % (2**n)
        return from_int(next_val)

    def D_op(S):
        """The decomposition operator D."""
        res = set()
        for s in S:
            res.update(s)
        return res

    def C_op(sigma):
        """The composition operator C."""
        # For each component k, find the set of possible values from sigma.
        # This is the core of the C operator.
        component_options = []
        for k in range(n):
            options_k = sigma.intersection(V[k])
            if not options_k:
                # Rule 1: if no values are specified, all are possible.
                component_options.append(list(V[k]))
            else:
                component_options.append(list(options_k))
        
        # The result is the Cartesian product of the component options.
        return set(itertools.product(*component_options))

    # ---- Simulation ----
    print(f"Running simulation for an {n}-bit counter.")
    print("This will show the size of the set C(sigma), which dictates the number of f() calls per step.")

    # Initial state: 0
    s0 = from_int(0)
    
    # Initial sigma set for relaxed simulation
    sigma = D_op({s0})
    
    # An ordinary simulation would just compute f(s_i) once per step.
    # Cost of ordinary simulation per step = 1 f() call.

    for i in range(n):
        # Calculate C(sigma_i)
        C_sigma = C_op(sigma)
        
        # The size of this set is the number of f() calls we must make.
        cost = len(C_sigma)
        
        print(f"Step {i}:")
        print(f"  - Size of value set |sigma_{i}|: {len(sigma)}")
        print(f"  - Cost of relaxed step: |C(sigma_{i})| = {cost} evaluations of f")
        
        if i < n -1: # Avoid huge computation on last step
            # Calculate the next sigma set
            next_states_D = set()
            for s_intermediate in C_sigma:
                s_next = f(s_intermediate)
                next_states_D.update(D_op({s_next}))
            
            sigma.update(next_states_D)
    
    print("\nConclusion:")
    print("The cost of the relaxed simulation (number of f() evaluations) grows as 2, 4, 8, ...")
    print("This is exponential growth compared to the ordinary simulation's constant cost of 1.")
    print("This demonstrates that for certain functions, the relaxed simulation requires")
    print("exponentially larger computational effort and memory to hold the intermediate states.")


demonstrate_exponential_blowup()