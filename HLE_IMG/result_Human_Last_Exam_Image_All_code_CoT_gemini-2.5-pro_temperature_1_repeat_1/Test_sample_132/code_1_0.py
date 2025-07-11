def calculate_conversation_probability():
    """
    This function calculates the probability of a specific conversation path
    in a fixed Markov chain model.
    """
    # Probabilities of the transitions in the identified path
    # Path: START -> Greeting -> Promise -> Solution -> END
    
    # P(START -> Greeting)
    p_start_greeting = 0.28
    
    # P(Greeting -> Promise)
    p_greeting_promise = 0.53
    
    # P(Promise -> Solution), which is the missing probability from the 'Promise' state (1.0 - 0.74 - 0.21)
    p_promise_solution = 0.05
    
    # P(Solution -> END), as given in the problem description
    p_solution_end = 1.0
    
    # The total probability of the path is the product of individual transition probabilities
    total_probability = p_start_greeting * p_greeting_promise * p_promise_solution * p_solution_end
    
    # Print the final equation with all the numbers
    print("The conversation maps to the path: START -> Greeting -> Promise -> Solution -> END")
    print("The probability is the product of the transition probabilities in this path.")
    print(f"P = {p_start_greeting} * {p_greeting_promise} * {p_promise_solution} * {p_solution_end}")
    
    # Print the final rounded result
    print(f"Calculated Probability = {total_probability:.4f}")

calculate_conversation_probability()