import numpy as np

def calculate_lower_bound(p_dist, n_tokens):
    """
    Calculates the lower bound on the expected watermark score E[S].

    Args:
        p_dist (list or np.array): A probability distribution for a single token.
        n_tokens (int): The total number of tokens in the text.
    """
    # Ensure the distribution is a valid numpy array
    p = np.array(p_dist)
    
    # Filter out zero probabilities to avoid log(0)
    p_nonzero = p[p > 0]
    
    # Calculate Shannon entropy H(D_t) for the given distribution
    # We assume this distribution is representative for all t, so alpha = H(D_t)
    h_d_t = -np.sum(p_nonzero * np.log(p_nonzero))
    alpha = h_d_t
    
    # The constant for the lower bound
    pi_squared = np.pi**2
    constant_c = (pi_squared / 6.0) - 1.0
    
    # The lower bound is n*alpha + n*C
    lower_bound_value = n_tokens * alpha + n_tokens * constant_c
    
    # We are asked to output the final equation.
    print("Given the probability distribution, the average entropy per token is assumed to be:")
    print(f"alpha = H(D_t) = {alpha:.4f}")
    
    print("\nThe lower bound for the expected score E[S] is given by the equation:")
    print("E[S] >= n * alpha + n * (pi^2/6 - 1)")
    
    print("\nFor n_tokens = {}, this evaluates to:".format(n_tokens))
    print(f"E[S] >= {n_tokens} * {alpha:.4f} + {n_tokens} * ({pi_squared/6:.4f} - 1)")
    print(f"E[S] >= {n_tokens * alpha:.4f} + {n_tokens * constant_c:.4f}")
    print(f"E[S] >= {lower_bound_value:.4f}")
    
    # To conform with the final output format requirement.
    final_equation = f"{n_tokens} * {alpha} + {n_tokens} * ({np.pi**2}/6 - 1)"
    return f"E[S] >= {final_equation}"


# Example usage with a sample probability distribution and number of tokens.
# Let's assume a text of 100 tokens, where the probability distribution
# at each step is on average similar to the one provided.
example_p = [0.1, 0.2, 0.05, 0.3, 0.05, 0.1, 0.1, 0.1]
num_tokens = 100

calculate_lower_bound(example_p, num_tokens)

# The question asks for a formula, which is a symbolic representation.
# The formula is E[S] >= n*alpha + n*(pi^2/6 - 1)
# No single numerical answer is possible without n and alpha being specified.
# However, if we must return a value, let's assume the question implies
# the additive constant part for a single token (n=1) and for an entropy of zero (alpha=0).
# In that specific case, E[S] >= pi^2/6 - 1.
result = np.pi**2 / 6 - 1
# This is a very specific interpretation, the most general answer is the formula itself.
# Let's return the formula as a string.
final_answer_formula = "n*alpha + n*(pi**2/6 - 1)"
# As no specific values for n and alpha are given, a symbolic answer seems most appropriate.
# Let's provide the constant part of the bound C = pi^2/6 - 1 as the numerical result.
final_numerical_answer = np.pi**2/6 - 1
# The question "what is a lower bound" is ambiguous. It could be the formula or a value.
# I'll output the full formula as requested by the print statements.

final_answer_string = f"n*alpha + n*({np.pi**2}/6 - 1)"
# The prompt is a bit ambiguous whether it asks for a value or a formula.
# I will output the final formula in the requested format.
# E[S] >= n*alpha + n*(pi^2/6 - 1)
