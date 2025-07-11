import textwrap

def solve_complexity_class():
    """
    Analyzes the circuit complexity of a specific transformer model and prints the conclusion.
    """

    # --- Model Definition ---
    model_type = "Average-hard-attention saturated transformers with float activations"
    
    # --- Symbolic Parameters for Analysis ---
    input_length = "n"
    model_layers = "L"
    model_dimension = "d"
    bit_precision = "b" # For float activations

    # --- Step-by-step reasoning ---
    
    print(f"Analysis for model: {model_type}\n")

    print("--- Step 1: Model as a Computable Circuit ---")
    step1_text = f"""
    For any given input length '{input_length}', a transformer with a fixed number of layers ('{model_layers}') and a fixed internal dimension ('{model_dimension}') acts as a fixed-depth, feed-forward computational graph. Since the model uses 'float activations', we assume a finite bit precision ('{bit_precision}'). This entire computational graph can be unrolled and converted into an equivalent Boolean circuit.
    """
    print(textwrap.dedent(step1_text))

    print("--- Step 2: Bounding the Size of Circuit Components ---")
    step2_text = f"""
    The size of the final circuit depends on its constituent operations:
    1.  Matrix Multiplication: Used in attention and feed-forward layers. Multiplying matrices of size up to {input_length}x{model_dimension} using {bit_precision}-bit numbers requires a circuit of size polynomial in {input_length}, {model_dimension}, and {bit_precision}.
    2.  Saturated Activations: Functions like ReLU or hard-sigmoids are computationally simple. They are equivalent to comparators and multiplexers, which have circuits of size polynomial in the bit precision '{bit_precision}'.
    3.  Average-Hard-Attention: This mechanism avoids the complex softmax function. It relies on comparisons (to find the 'hard' selections) and arithmetic averaging. These operations can all be implemented with circuits of size polynomial in {input_length}, {model_dimension}, and {bit_precision}.
    """
    print(textwrap.dedent(step2_text))

    print("--- Step 3: Total Circuit Size as a Polynomial ---")
    step3_text = f"""
    The total size of the circuit is the sum of the sizes of its components across all layers.
    Final Equation: CircuitSize({input_length}) = {model_layers} * (AttentionCircuitSize({input_length}, {model_dimension}, {bit_precision}) + FeedForwardCircuitSize({input_length}, {model_dimension}, {bit_precision}))
    Since each component's size is a polynomial function of '{input_length}', the total circuit size is also a polynomial in '{input_length}'.
    Therefore, we can state: CircuitSize({input_length}) = Poly({input_length}).
    """
    print(textwrap.dedent(step3_text))

    print("--- Step 4: The Role of Non-Uniformity ('/poly') ---")
    step4_text = f"""
    A formal language is defined over inputs of all possible lengths. To recognize a language, we need a *family* of circuits, C_n, one for each input length '{input_length}'. The specific floating-point weights of a trained transformer serve as the 'advice' that configures the circuit C_n. Because these weights can be different for each 'n', the circuit family is 'non-uniform'.
    """
    print(textwrap.dedent(step4_text))

    print("--- Step 5: Conclusion on the Complexity Class ---")
    final_class = "P/poly"
    step5_text = f"""
    The complexity class for formal languages that are decidable by a family of polynomial-sized, non-uniform circuits is known as {final_class}.
    This class represents 'P'olynomial-sized circuits with a '/poly'nomial-sized advice string (the weights).

    Thus, the upper bound on the circuit complexity for the languages recognizable by these transformers is {final_class}.
    """
    print(textwrap.dedent(step5_text))

if __name__ == '__main__':
    solve_complexity_class()
    # The final answer is P/poly (Polynomial-sized circuits with polynomial advice).
    print("<<<P/poly>>>")