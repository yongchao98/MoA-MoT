import math

def calculate_probability(n_vars, poly_degree, closeness_factor):
    """
    Calculates an upper bound on the probability that a random PCP proof is
    'close' to the set of correct proofs for a YES-instance.

    Args:
        n_vars (int): The number of variables in the original NP problem (e.g., 3-SAT).
        poly_degree (int): The proof length is n_vars^poly_degree.
        closeness_factor (float): A factor < 0.5 defining 'close'.
    """
    # The length of the proof string in a PCP system.
    # It's polynomial in the original problem size 'n_vars'.
    m_proof_len = n_vars ** poly_degree

    # The number of witnesses (e.g., satisfying assignments) is at most 2^n_vars.
    # This gives an upper bound on the number of correct proofs |Π(x)|.
    k_max_proofs = 2**n_vars

    # We want to calculate the probability that a random proof `π_rand` has a relative
    # distance less than `closeness_factor` to a specific correct proof `π*`.
    # Let d = closeness_factor * m_proof_len.
    # P(δ(π_rand, π*) < closeness_factor) = (1/2^m) * Σ_{i=0 to d} C(m, i)
    # where C(m, i) is "m choose i".
    d_max_dist = int(closeness_factor * m_proof_len)

    try:
        # Sum of binomial coefficients C(m, i) for i from 0 to d.
        sum_of_combs = sum(math.comb(m_proof_len, i) for i in range(d_max_dist + 1))
    except (ValueError, OverflowError) as e:
        print(f"Could not compute for n={n_vars}, d={poly_degree}: {e}")
        return

    # Probability of being close to a *single* correct proof.
    p_single_proof = sum_of_combs / (2**m_proof_len)

    # Use the union bound to find the upper bound on the probability of being close
    # to *any* proof in the set Π(x).
    # P(close to Π(x)) <= |Π(x)| * P(close to a single π*)
    p_total_bound = k_max_proofs * p_single_proof

    print("--- Illustrating the BPP Algorithm for NP ---")
    print(f"Problem variables n = {n_vars}")
    print(f"Proof length m = n^{poly_degree} = {m_proof_len}")
    print(f"Max number of correct proofs |Π(x)| <= 2^{n_vars} = {k_max_proofs}")
    print(f"Defining 'close' as relative distance < {closeness_factor}\n")

    print("The probability that a random proof is 'close' to a single specific correct proof is:")
    # Using scientific notation for small numbers.
    print(f"P(δ(π_rand, π*) < {closeness_factor}) = {p_single_proof:.2e}\n")

    print("The final upper bound on the probability that a random proof is close to *any* correct proof is:")
    equation_str = (f"P(δ(π_rand, Π(x)) < {closeness_factor}) "
                    f"<= |Π(x)| * P(δ(π_rand, π*) < {closeness_factor})")
    print(equation_str)

    # We output each number in the final equation as requested
    final_calc_str = (f"                                "
                      f"<= {k_max_proofs} * {p_single_proof:.2e} "
                      f"= {p_total_bound:.2e}")
    print(final_calc_str)
    print("\nThis probability is negligible, so we can reliably distinguish YES instances (where δ ≈ 0.5) from NO instances (where δ = 1).")


# --- Parameters for the demonstration ---
# A small problem size
problem_size_n = 20
# A conservative polynomial degree for proof length (e.g., m = n^2)
proof_length_poly_degree = 2
# A generous definition of 'close' (anything less than 40% distance)
distance_factor = 0.4

calculate_probability(problem_size_n, proof_length_poly_degree, distance_factor)