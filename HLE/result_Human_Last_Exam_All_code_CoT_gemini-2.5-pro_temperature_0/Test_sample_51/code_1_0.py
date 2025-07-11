def demonstrate_paradox():
    """
    This script illustrates the logical inconsistency that arises from combining
    general recursion with the Law of Excluded Middle (LEM).

    The problem's subterm rule allows for general recursion, meaning we can write
    any computable function. LEM (Axiom H) asserts P ∨ ¬P for any proposition P.
    Applying LEM to the Halting Problem implies the existence of a function
    `decide_halts`, which is known to be impossible to create. This script shows
    the contradiction that arises from assuming such a function exists.
    """
    print("--- The Halting Problem Paradox ---")
    print("We assume a function 'decide_halts(p_source)' exists, as implied by LEM.")
    print("This function predicts if a program `p` halts on its own source `p_source`.\n")

    def paradox_program_factory(decide_halts_func):
        """
        This factory creates a paradoxical program based on a given halt decider.
        The program does the opposite of what the decider predicts.
        """
        def paradox(p_source_code):
            print(f"--- Running paradox on program source: '{p_source_code}' ---")
            # Use the decider to predict the behavior of p(p)
            prediction = decide_halts_func(p_source_code)
            print(f"Halt decider predicts that '{p_source_code}' will halt: {prediction}")

            if prediction is True:
                # If it's predicted to halt, we do the opposite and loop forever.
                print("Prediction is HALT, so I will now loop forever.")
                print("Actual behavior: DOES NOT HALT")
                return "DOES NOT HALT"
            else:
                # If it's predicted to not halt, we do the opposite and halt immediately.
                print("Prediction is DO_NOT_HALT, so I will now halt.")
                print("Actual behavior: HALTS")
                return "HALTS"
        return paradox

    paradox_source = "paradox_program"

    print("Analyzing the call: paradox(paradox_program)")
    print("We examine the two cases for the decider's prediction.\n")

    # Case 1: Assume the halt decider predicts that paradox(paradox_program) will halt.
    print("="*60)
    print("CASE 1: Assume decide_halts('paradox_program') returns True (predicts HALT).")
    print("="*60)

    def decider_predicts_halt(source):
        if source == paradox_source:
            return True
        return False # Default for other programs

    paradox_v1 = paradox_program_factory(decider_predicts_halt)
    result1 = paradox_v1(paradox_source)
    print("\nCONTRADICTION: The program was predicted to HALT, but its actual behavior is DOES NOT HALT.")
    print("="*60)
    print("\n")

    # Case 2: Assume the halt decider predicts that paradox(paradox_program) will not halt.
    print("="*60)
    print("CASE 2: Assume decide_halts('paradox_program') returns False (predicts NOT HALT).")
    print("="*60)

    def decider_predicts_no_halt(source):
        if source == paradox_source:
            return False
        return True # Default for other programs

    paradox_v2 = paradox_program_factory(decider_predicts_no_halt)
    result2 = paradox_v2(paradox_source)
    print("\nCONTRADICTION: The program was predicted to NOT HALT, but its actual behavior is HALTS.")
    print("="*60)
    print("\n")

    print("Conclusion: The existence of a `decide_halts` function leads to a logical contradiction.")
    print("Since the Law of Excluded Middle (H) implies this function's existence in a system with general recursion, the combination is inconsistent.")

if __name__ == '__main__':
    demonstrate_paradox()