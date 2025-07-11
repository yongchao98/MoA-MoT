import textwrap

def find_jzdec_formula():
    """
    This function explains the derivation of the Linear Logic formula for a
    Minsky machine's JZDEC instruction and prints the result.
    """

    # Introduction to the components of the formula.
    # S_l: A literal representing that the machine is in state l.
    # C_r: A literal representing a single unit in counter r.
    # Z_r: A special literal for verifying that counter r is zero.
    # r, z, nz: Symbolic placeholders for the instruction's parameters.
    #           r = register index, z = zero-branch state, nz = non-zero-branch state.

    # The logical connectives used:
    # --o : Linear Implication (multimap)
    # (x) : Multiplicative Conjunction, Tensor Product (otimes)
    # &   : Additive Conjunction, 'with' (internal choice)

    # Step 1: Model the non-zero case (x_r > 0).
    # The instruction decrements counter r and jumps to state nz.
    # In Linear Logic, this means we must consume one C_r resource and
    # produce one S_nz resource. The linear implication `C_r --o S_nz`
    # perfectly models this: given a C_r, it yields an S_nz.
    # This path is only provable if a C_r resource is available.
    
    # We use strings to represent the symbolic formula components.
    r_sym, z_sym, nz_sym = 'r', 'z', 'nz'
    
    f_nonzero = f"(C_{r_sym} --o S_{nz_sym})"

    # Step 2: Model the zero case (x_r = 0).
    # The instruction jumps to state z. This must only be possible if x_r = 0.
    # The problem provides a mechanism for this: the Z_r literal. The axioms in Delta
    # are designed so that Z_r can only be fully resolved (i.e., eliminated
    # from the proof context) if no C_r resources are present.
    # Therefore, the formula for this case must produce the new state S_z and
    # also produce the Z_r literal to trigger the zero-test.
    # Producing both resources at once is modeled by the tensor product.
    
    f_zero = f"(S_{z_sym} (x) Z_{r_sym})"

    # Step 3: Combine the two cases.
    # The final formula F(r, z, nz) must offer a choice between the two cases.
    # When proving the sequence, the prover will be forced to choose the branch
    # that is consistent with the available resources (the counter state).
    # For example, if C_r is present, the prover can choose the f_nonzero branch.
    # If C_r is absent, that branch fails, and the prover must choose the f_zero branch.
    # This type of choice is modeled by the '&' (with) connective.
    
    final_formula = f"{f_nonzero} & {f_zero}"

    # Print the final result and explanation.
    # The prompt asks to "output each number in the final equation".
    # Since the formula is symbolic, we print the symbolic names 'r', 'z', and 'nz'.
    
    explanation = textwrap.dedent(f"""
    The appropriate formula F(r, z, nz) that encodes the JZDEC instruction is derived
    by combining the logic for its two possible outcomes.

    1.  Non-Zero Case (x_r > 0): Decrement counter r and go to state nz.
        This is encoded as: {f_nonzero}
        This consumes a C_r resource and produces an S_nz resource.

    2.  Zero Case (x_r = 0): Go to state z.
        This is encoded as: {f_zero}
        This produces an S_z resource and triggers the zero-test for counter r with Z_r.

    These two mutually exclusive outcomes are presented as an internal choice using the
    '&' (with) connective. The correct path is determined by the context (the value of counter r).

    Final Formula F({r_sym}, {z_sym}, {nz_sym}):
    """).strip()

    print(explanation)
    print(final_formula)
    
    legend = textwrap.dedent("""
    
    Notation:
      --o : Linear Implication (multimap)
      (x) : Tensor Product (otimes)
      &   : 'With' connective (internal choice)
    """).strip()
    
    print(legend)

find_jzdec_formula()