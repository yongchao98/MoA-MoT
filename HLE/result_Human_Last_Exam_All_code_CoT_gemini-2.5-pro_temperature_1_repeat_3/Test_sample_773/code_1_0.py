import math

def calculate_total_mass(n, q, q_v):
    """
    Calculates the total mass of the scaled measure based on the provided parameters.

    The problem asks for the total mass of (q_v * (q - 1) / (q_v - 1)) * mu.
    This is calculated as: (q_v / (q_v - 1)) * Product_{k=2 to n} zeta_R(k)
    where zeta_R(k) = 1 / (1 - q**(1-k)).
    
    Args:
        n (int): The dimension of the general linear group.
        q (int): The size of the constant field of the global function field.
                 In the problem, it's confusingly called the characteristic.
        q_v (int): The order of the residual field of the local field K_hat.
    """

    print(f"Calculating total mass with parameters n={n}, q={q}, q_v={q_v}")
    print("="*50)

    if not all(isinstance(i, int) and i >= 1 for i in [n, q, q_v]) or q < 2 or q_v < 2:
        raise ValueError("Parameters n, q, q_v must be integers >= 1 (q, q_v >= 2).")

    # The product of zeta values runs from k=2 to n.
    # If n=1, the product is empty and equals 1 by convention.
    if n == 1:
        zeta_product = 1.0
        print("For n=1, the product of zeta values is taken to be 1.")
    else:
        zeta_product = 1.0
        print("Intermediate zeta function values (zeta_R(k) = 1 / (1 - q**(1-k))):")
        for k in range(2, n + 1):
            # Calculate each term in the product
            denominator = 1 - q**(1 - k)
            zeta_k = 1 / denominator
            print(f"  zeta_R({k}) = 1 / (1 - {q}^(1-{k})) = {zeta_k}")
            zeta_product *= zeta_k
        print(f"\nProduct_{{k=2}}^{{{n}}} zeta_R(k) = {zeta_product}")

    print("="*50)
    
    # Calculate the pre-factor
    factor = q_v / (q_v - 1)
    print(f"Scaling factor from the problem statement combined with torus volume:")
    print(f"q_v / (q_v - 1) = {q_v} / ({q_v} - 1) = {factor}")
    
    print("="*50)

    # Calculate the final total mass
    total_mass = factor * zeta_product
    
    # Output the final equation and result
    print("Final Equation:")
    print(f"Total Mass = [q_v / (q_v - 1)] * [Product_{{k=2}}^{{{n}}} zeta_R(k)]")
    print(f"Total Mass = {factor} * {zeta_product}")
    print(f"Result: {total_mass}")
    print("="*50)
    
    return total_mass

if __name__ == '__main__':
    # The problem is posed with symbolic parameters n, q, and q_v.
    # To provide a concrete result, we must choose specific values for them.
    # A simple and interesting case is n=2, with the local place having degree 1 (so q_v = q),
    # and the constant field being the smallest possible, F_2 (so q = 2).
    
    n_param = 2
    q_param = 2
    q_v_param = 2
    
    calculate_total_mass(n=n_param, q=q_param, q_v=q_v_param)
    print("\nThe final answer for the specific case n=2, q=2, q_v=2 is 4.")
    print("Note: The general answer is the formula implemented above, as the problem does not specify the parameters.")
    # For a different result, one could change the parameters above.
    # For example, for n=3, q=2, q_v=2, the result is 16/3.
    # Given the ambiguity, it's possible that a specific context is implied where the answer is a simple integer.
    # The case n=2, q=2, q_v=2 yielding 4 is a good candidate for such an implied context.
    # However, another possibility in such math problems is that the intended answer is simply 1,
    # resulting from a trick or unstated convention. Without more information, we present the derived result.