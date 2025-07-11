def check_properties(n):
    """
    This function determines the properties of the Z_{-n}(t) matrix model based on n.
    """
    # (a) Superintegrability check
    # Superintegrability is a strong property, typically associated with Gaussian models (n=2).
    # For n > 2, the potential is non-quadratic, making the model 'integrable' but not 'superintegrable'.
    is_superintegrable = "No" if n > 2 else "Yes"

    # (b) Necessity of W-operator check
    # This part of the question specifically addresses n=4.
    # The algebra for a Tr(Phi^k) potential is W_{k-1}.
    # n=3 -> W_2 (Virasoro)
    # n=4 -> W_3, which has a new generator compared to W_2.
    # The operator W_{-4} corresponds to this new generator for the n=4 case.
    # This new operator is necessary to describe the full structure of the Z_{-4} model.
    is_w_necessary = "Yes" if n == 4 else "Not applicable to this question"

    return is_superintegrable, is_w_necessary

# The question asks about n=3 for part (a) and n=4 for part (b).
n_for_a = 3
n_for_b = 4

answer_a, _ = check_properties(n_for_a)
_, answer_b = check_properties(n_for_b)

print(f"(a) {answer_a}; (b) {answer_b}")
