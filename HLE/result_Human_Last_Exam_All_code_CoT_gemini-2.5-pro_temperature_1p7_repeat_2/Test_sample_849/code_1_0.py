import sys

# A map of Gödel numbers (indices) to predictor names for clarity
predictor_names = {
    1: "Alternator (predicts 1 - last_digit)",
    2: "Repeater (predicts last_digit)",
    3: "Fixed Pattern '110' (predicts based on position)",
    4: "Always Zero (always predicts 0)"
}

# Define the predictors (our partial computable functions)
def predictor_alternator(sequence):
    """Predicts the opposite of the last digit."""
    if not sequence:
        return 0 # Default prediction for empty sequence
    return 1 - sequence[-1]

def predictor_repeater(sequence):
    """Predicts the same as the last digit."""
    if not sequence:
        return 1 # Default prediction for empty sequence
    return sequence[-1]

def predictor_fixed_pattern(sequence):
    """Predicts based on a fixed pattern 1,1,0,1,1,0,..."""
    pattern = [1, 1, 0]
    return pattern[len(sequence) % 3]

def predictor_always_zero(sequence):
    """Always predicts 0, a very simple hypothesis."""
    return 0

# Store the predictors in a dictionary mapping index to function
predictors = {
    1: predictor_alternator,
    2: predictor_repeater,
    3: predictor_fixed_pattern,
    4: predictor_always_zero
}

class ModelM:
    """
    A model that pursues a safe goal by performing verifiable learning from
    logical deductions. It identifies predictors that are consistent with
    the known sequence, rather than guessing which one will succeed.
    """
    def __init__(self, predictor_map):
        self.predictors = predictor_map
        self.predictor_indices = list(predictor_map.keys())

    def _is_consistent(self, predictor_func, sequence):
        """
        Checks if a predictor's outputs are logically consistent
        with the observed sequence. This is a verifiable, deductive step.
        A predictor is consistent if it correctly "predicts" each digit
        of the sequence given the part of the sequence that came before it.
        """
        # We check from t=1 to the end of the sequence.
        for t in range(len(sequence)):
            # The context is the sequence up to time t
            context = sequence[:t]
            # The actual digit at time t
            actual_digit = sequence[t]
            # The predictor's output for that context
            predicted_digit = predictor_func(context)
            # If the prediction is wrong, the predictor is inconsistent
            if predicted_digit != actual_digit:
                return False
        return True

    def find_consistent_predictors(self, initial_sequence):
        """
        The main safe goal: Instead of trying to guess the correct predictor for
        the future, this method identifies the set of all predictors that
        have not yet been falsified by the evidence.
        """
        consistent_set = set()
        print(f"Analyzing initial sequence: {initial_sequence}")
        print("Verifying consistency of all predictors...")
        print("-" * 40)
        for index in sorted(self.predictor_indices):
            p_func = self.predictors[index]
            p_name = predictor_names[index]
            if self._is_consistent(p_func, initial_sequence):
                consistent_set.add(index)
                print(f"  [PASS] Predictor {index} ({p_name}): CONSISTENT")
            else:
                print(f"  [FAIL] Predictor {index} ({p_name}): INCONSISTENT")
        print("-" * 40)
        return consistent_set

def main():
    # The user provides an initial part of a binary sequence
    # Example: A sequence that starts as an alternating series
    initial_sequence = [0, 1, 0, 1]

    # The model M is instantiated with its hypothesis classes (predictors)
    model = ModelM(predictors)

    # M's safe goal is to find the set I_c of predictors consistent with the data.
    # This is the result of "Verifiable learning from logical deductions".
    consistent_indices = model.find_consistent_predictors(initial_sequence)

    # Output the final result based on the safe goal.
    # This represents the model's safe state of knowledge.
    # It hasn't proven which predictor is correct for the future, but it has
    # proven which ones are logically consistent with the available evidence.
    print("\nThe safe goal's output is not a single prediction, but the set of")
    print("all hypotheses that are logically consistent with the known data.\n")
    
    print("Final Equation (Set of Consistent Predictor Indices):")
    # Print the equation part by part to fulfill the "output each number" requirement
    equation_lhs = "I_consistent = {"
    sys.stdout.write(equation_lhs)
    
    # Iterate and print each number in the set
    num_list = sorted(list(consistent_indices))
    for i, num in enumerate(num_list):
        sys.stdout.write(str(num))
        if i < len(num_list) - 1:
            sys.stdout.write(", ")
            
    equation_rhs = "}"
    sys.stdout.write(equation_rhs + "\n")


if __name__ == '__main__':
    main()
