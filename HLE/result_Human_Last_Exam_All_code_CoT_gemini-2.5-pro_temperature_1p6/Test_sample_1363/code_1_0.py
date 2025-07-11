import textwrap

def solve_dance_question():
    """
    Analyzes the technical principles of ballroom dances to determine where
    overturning a reverse turn is fundamentally against the technique.
    """
    options = {
        'A': 'Viennese Waltz',
        'B': 'English Waltz',
        'C': 'European Tango',
        'D': 'Slow Foxtrot',
        'E': 'Quickstep'
    }

    print("### Analyzing the Possibility of Overturning a Reverse Turn ###")
    print("-" * 60)
    
    analysis_swing = """
    In the 'swing' dances (Waltz, Viennese Waltz, Foxtrot, Quickstep), the technique is based on principles of CBM (Contrary Body Movement), body swing, and sway. These elements provide a great deal of rotational freedom. While figures have a standard amount of turn, it is common and technically permissible in advanced choreography to 'overturn' or 'underturn' them to create different effects or to navigate the dance floor. This flexibility is part of the dance's character.
    """
    
    analysis_tango = """
    The European Tango is fundamentally different. Its technique has no body swing, no sway, and no rise and fall. The movement is staccato, sharp, and grounded. A critical technical aspect is 'No Foot Turn' on many steps, where rotation is generated from the body over the feet. To overturn a Reverse Turn in Tango, a dancer would be forced to introduce sway or pivot on a foot in a way that is not allowed, thereby breaking the core, established technique of the dance.
    """

    print("1. Analysis of Swing Dances (A, B, D, E):")
    print(textwrap.fill(analysis_swing, width=70))
    print("\n2. Analysis of European Tango (C):")
    print(textwrap.fill(analysis_tango, width=70))
    print("-" * 60)

    print("### Conclusion presented as a logical equation: ###")
    print("Let 'is_overturn_possible' be a function that returns True if a reverse turn can be overturned without breaking technique.")
    print("")
    print("is_overturn_possible(A. Viennese Waltz) = True")
    print("is_overturn_possible(B. English Waltz)  = True")
    print("is_overturn_possible(D. Slow Foxtrot)   = True")
    print("is_overturn_possible(E. Quickstep)      = True")
    print("is_overturn_possible(C. European Tango) = False")
    print("")

    answer_key = 'C'
    print(f"The only dance for which the statement is False is {options[answer_key]}.")

solve_dance_question()