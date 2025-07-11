import numpy as np

def calculate_entropy(state_distribution):
    """
    Calculates the entropy H(p) = - sum(p * log2(p)) for a given state distribution.
    This is the objective function that the policy optimizes in the limit k->infinity.
    """
    p = np.array(state_distribution)
    # Add a small epsilon to avoid log(0) for states with zero probability
    entropy_val = -np.sum(p * np.log2(p + 1e-12))
    return entropy_val

def final_objective_equation_demo():
    """
    Demonstrates the calculation of the objective function at the limit k->infinity,
    which is the state entropy H(s).
    """
    print("The iterative process converges to a policy that maximizes its own state entropy.")
    print("The objective at the limit is: J(pi_inf) = max H(s) = max [-sum(p(s) * log(p(s)))]")

    # For a system with a finite number of states, maximum entropy is achieved
    # when the state distribution is uniform.
    num_states = 4
    # The maximum entropy policy results in a uniform state distribution.
    uniform_distribution = [1/num_states] * num_states
    
    print(f"\nFor a system with {num_states} states, the maximum entropy is achieved with a uniform distribution:")
    print(f"p(s) = {uniform_distribution}")
    
    entropy_value = calculate_entropy(uniform_distribution)

    print("\nThe final equation being maximized for this uniform distribution is:")
    
    # We build and print the equation string with each term.
    equation_str = "H(s) = - ["
    terms = []
    for p_s in uniform_distribution:
        # Each term in the sum is p(s) * log2(p(s))
        term_val = p_s * np.log2(p_s)
        terms.append(f"({p_s:.2f} * log2({p_s:.2f}))")
    equation_str += " + ".join(terms) + "]"
    print(equation_str)
    
    print("\nCalculating each term:")
    p_s = uniform_distribution[0]
    log_p_s = np.log2(p_s)
    term_val = p_s * log_p_s
    print(f"Each term is {p_s:.2f} * {log_p_s:.2f} = {term_val:.2f}")
    print(f"Total H(s) = - [4 * {term_val:.2f}] = {entropy_value:.4f} bits")
    
    print("\nThis maximum entropy value is the objective achieved by the limiting policy lim_{k->inf} pi^k.")

final_objective_equation_demo()