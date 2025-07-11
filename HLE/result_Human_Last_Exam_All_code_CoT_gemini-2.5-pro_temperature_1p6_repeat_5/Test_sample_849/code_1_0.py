import sys

# Step 1: Define the ground truth sequence (Thue-Morse sequence)
# The rule is: T(0) = 0; T(2n) = T(n); T(2n+1) = 1 - T(n).
def generate_thue_morse(n):
    """Generates the Thue-Morse sequence of length n."""
    sequence = [0]
    while len(sequence) < n:
        sequence.extend([1 - x for x in sequence])
    return sequence[:n]

# Step 2: Define the predictors (p.c. functions)
# Each predictor takes a sequence prefix and returns a prediction (0 or 1).

def predictor_zero(prefix):
    """Always predicts 0."""
    return 0

def predictor_one(prefix):
    """Always predicts 1."""
    return 1

def predictor_alternate(prefix):
    """Predicts the opposite of the last digit."""
    if not prefix: return 0
    return 1 - prefix[-1]

def predictor_repeat(prefix):
    """Predicts the same as the last digit."""
    if not prefix: return 0
    return prefix[-1]
    
def predictor_correct_thue_morse(prefix):
    """This predictor 'knows' the Thue-Morse rule."""
    n = len(prefix)
    if n == 0:
        return 0
    # Use the recursive definition to find the next element
    if n % 2 == 0:
        # T(2k) = T(k). The next element is at index n, so we need T(n).
        # We need the value at index n/2 from the *prefix*.
        return prefix[n // 2]
    else:
        # T(2k+1) = 1 - T(k). The next element is at index n, so we need 1 - T((n-1)/2).
        return 1 - prefix[(n - 1) // 2]


def main():
    """Main function to run the simulation."""
    # Simulation parameters
    sequence_length = 20
    initial_prefix_length = 3
    
    # Generate the full sequence
    ground_truth_sequence = generate_thue_morse(sequence_length)
    print(f"Ground Truth Sequence (Thue-Morse): {ground_truth_sequence}\n")

    # List of available predictors
    predictors = [
        predictor_zero,
        predictor_one,
        predictor_alternate,
        predictor_repeat,
        predictor_correct_thue_morse
    ]

    # The set I for indices of successful predictors
    # In this simulation, this set is computable because the sequence is finite.
    # In the general case, this could be uncomputable.
    index_set_I = set()

    # The initial part of the sequence known to the model
    initial_prefix = ground_truth_sequence[:initial_prefix_length]
    print(f"Starting with initial prefix: {initial_prefix}\n")
    print("-" * 50)

    # Step 3 & 4: Simulate prediction one step at a time
    # The model predicts the completion from index `initial_prefix_length` onwards.
    for t in range(initial_prefix_length, sequence_length):
        current_prefix = ground_truth_sequence[:t]
        actual_next_digit = ground_truth_sequence[t]
        
        print(f"Time t={t}:")
        print(f"  Prefix: {current_prefix}")
        print(f"  Actual next digit: {actual_next_digit}")
        print(f"  Predictions:")

        # For each predictor, get its prediction and check if it's correct.
        for i, predictor in enumerate(predictors):
            try:
                prediction = predictor(current_prefix)
                is_correct = (prediction == actual_next_digit)
                if is_correct:
                    # If correct, add the predictor's index to set I
                    index_set_I.add(i)
                    result_msg = "CORRECT (index added to I)"
                else:
                    result_msg = "WRONG"
                
                print(f"    - Predictor {i} ({predictor.__name__}): {prediction} -> {result_msg}")

            except Exception as e:
                # A p.c. function might not halt or might error.
                print(f"    - Predictor {i} ({predictor.__name__}): ERROR -> {e}")
        
        print("-" * 50)
        
    # Step 5: Output the final results
    # The prompt asks to "output each number in the final equation!".
    # We interpret this as printing the final computed set I.
    final_indices = sorted(list(index_set_I))
    
    print("\nSIMULATION COMPLETE\n")
    print("The final index set I contains indices of predictors that were correct at least once.")
    # Printing each number in the final set
    print(f"Final Index Set I = {{", end="")
    print(*final_indices, sep=", ", end="")
    print("}")


if __name__ == "__main__":
    main()
