def calculate_failure_probability():
    """
    This function demonstrates the failure of the proposed strategy for case (B).
    It calculates the success probability for a sequence crafted by an adversary.
    """
    # The set of boxes Alice focuses on for her random choice.
    C = set(range(1, 11))
    
    # The adversary crafts a sequence S where the difference set F_S is {1, 2}.
    F_S = {1, 2}
    
    # Calculate the number of boxes in C that are also in F_S.
    C_intersect_F_S = C.intersection(F_S)
    
    num_intersect = len(C_intersect_F_S)
    size_C = len(C)
    
    # The success probability is 1 - (|C intersect F_S| / |C|)
    success_prob_numerator = size_C - num_intersect
    success_prob_denominator = size_C
    success_probability = success_prob_numerator / success_prob_denominator
    
    print("Let C be the set of boxes Alice randomly chooses from: {1, 2, ..., 10}")
    print(f"Let the adversary choose a sequence S where its difference set F_S = {F_S}")
    print(f"The set of failing guesses for Alice is C intersect F_S = {C_intersect_F_S}")
    print(f"The number of failing guesses is |C intersect F_S| = {num_intersect}")
    print(f"The size of C is |C| = {size_C}")
    print("Alice's success probability P is calculated as: P = (|C| - |C intersect F_S|) / |C|")
    print(f"P = ({size_C} - {num_intersect}) / {size_C} = {success_prob_numerator}/{success_prob_denominator} = {success_probability}")
    
    required_prob = 9/10
    print(f"The required success probability is at least {required_prob}")
    if success_probability < required_prob:
        print(f"Since {success_probability} < {required_prob}, the strategy is not guaranteed to succeed.")

calculate_failure_probability()
