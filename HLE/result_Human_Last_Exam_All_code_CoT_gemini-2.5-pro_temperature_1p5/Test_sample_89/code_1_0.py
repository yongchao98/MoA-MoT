import numpy as np

def calculate_learning_difficulty(p_values):
    """
    Calculates the information content and relative number of triplets for a given alignment p.

    Args:
        p_values (list or np.array): A list of alignment probabilities (p) between 0 and 1.
    """
    print("p: Probabilistic Representational Alignment")
    print("I(p): Information per triplet (in bits)")
    print("N: Relative number of triplets required (proportional to 1/I(p))")
    print("-" * 65)

    for p in p_values:
        # Edge cases for log(0)
        if p == 0 or p == 1:
            # lim_{x->0} x*log2(x) = 0. So I(0)=I(1)=1.
            info = 1.0
        elif p == 0.5:
            # I(0.5) = 1 + 0.5*log2(0.5) + 0.5*log2(0.5) = 1 - 1 = 0.
            info = 0.0
        else:
            # Information I(p) = 1 + p*log2(p) + (1-p)*log2(1-p)
            info = 1.0 + p * np.log2(p) + (1 - p) * np.log2(1 - p)

        # Number of triplets N is proportional to 1/I(p)
        if info < 1e-9:
            num_triplets = np.inf
        else:
            num_triplets = 1.0 / info
        
        # In the final equation N ~ 1/I, we output the numbers p, I, and N
        print(f"p = {p:<4.2f} | I(p) = {info:<6.4f} | N ~ 1/I(p) = {num_triplets:<8.4f}")


if __name__ == '__main__':
    # Demonstrate the relationship for key values of p
    demonstration_p_values = [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]
    calculate_learning_difficulty(demonstration_p_values)