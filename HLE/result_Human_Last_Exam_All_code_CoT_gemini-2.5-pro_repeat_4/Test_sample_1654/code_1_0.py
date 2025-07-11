def count_coat_descriptions():
    """
    Calculates the number of times Otchumyelov's coat is described
    to symbolize his shifting mentality in Chekov's "The Chameleon".
    """
    # In the story, Otchumyelov's interaction with his coat signals his changing mood.

    # He takes his coat off when he gets nervous and hot, fearing the dog's owner is important.
    times_coat_off = 2

    # He puts his coat on when he feels a "draught" and regains his composure,
    # believing the owner is a commoner.
    times_coat_on = 2
    
    # In the final scene, he wraps himself in his coat, representing his final disposition.
    final_wrap_action = 1

    # Calculate the total number of symbolic descriptions.
    total_descriptions = times_coat_off + times_coat_on + final_wrap_action

    print(f"Based on the analysis of the text, Chekov describes a change related to the coat {total_descriptions} times.")
    print(f"This is calculated from him taking it off {times_coat_off} times, putting it on {times_coat_on} times, and wrapping it around himself {final_wrap_action} time.")
    print(f"The final equation is: {times_coat_off} + {times_coat_on} + {final_wrap_action} = {total_descriptions}")

count_coat_descriptions()