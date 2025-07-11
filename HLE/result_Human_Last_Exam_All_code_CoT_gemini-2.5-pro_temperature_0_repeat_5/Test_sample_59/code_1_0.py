import numpy as np

def calculate_sbm_edge_probability():
    """
    Calculates the overall probability of an edge in a 2-block Stochastic Block Model (SBM).

    This serves as a concrete example for the general formula for jointly exchangeable graphs,
    which is P(edge) = integral of W(u,v) du dv, where W is the graphon.

    In this SBM example:
    - The graph has two communities (blocks) of equal size.
    - The latent space [0,1] is split at 0.5.
    - p_11 is the probability of an edge between two nodes in block 1.
    - p_22 is the probability of an edge between two nodes in block 2.
    - p_12 is the probability of an edge between a node in block 1 and a node in block 2.
    """
    # Define the parameters for our example SBM
    # Probability of connection within block 1 (nodes with latent vars in [0, 0.5))
    p_11 = 0.8
    # Probability of connection within block 2 (nodes with latent vars in [0.5, 1])
    p_22 = 0.7
    # Probability of connection between block 1 and block 2
    p_12 = 0.1

    # Proportions of the communities in the latent space [0,1]
    c1_size = 0.5
    c2_size = 1.0 - c1_size

    # The overall probability is the weighted average of the block probabilities,
    # where weights are the areas of the corresponding regions in the [0,1]x[0,1] square.
    # Area of (block 1, block 1) region = c1_size * c1_size
    # Area of (block 2, block 2) region = c2_size * c2_size
    # Area of (block 1, block 2) and (block 2, block 1) regions = 2 * c1_size * c2_size
    
    prob = (c1_size**2) * p_11 + (c2_size**2) * p_22 + 2 * c1_size * c2_size * p_12

    # Print the final equation with all the numbers, as requested.
    print(f"For a 2-block SBM with p_11={p_11}, p_22={p_22}, p_12={p_12} and block sizes {c1_size}, {c2_size}:")
    print(f"The probability to draw a link y_ij is P(y_ij=1)")
    print(f"= ({c1_size}**2) * {p_11} + ({c2_size}**2) * {p_22} + 2 * {c1_size} * {c2_size} * {p_12}")
    print(f"= {c1_size**2:.2f} * {p_11} + {c2_size**2:.2f} * {p_22} + {2*c1_size*c2_size:.2f} * {p_12}")
    print(f"= {c1_size**2 * p_11:.3f} + {c2_size**2 * p_22:.3f} + {2 * c1_size * c2_size * p_12:.3f}")
    print(f"= {prob:.3f}")

# Execute the function
calculate_sbm_edge_probability()