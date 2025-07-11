def solve_riddle():
    """
    This script identifies the thinker based on the provided context
    and prints the answer along with a calculation based on the years given.
    """
    
    # The years mentioned in the text for the visits.
    start_year = 1905
    end_year = 1909
    
    # The thinker in question.
    thinker_name = "Jacques Maritain"
    
    # The quote is from the memoirs of his wife, Raïssa, about their shared experiences.
    # The question asks for the "prominent figure of European thought," which best describes Jacques Maritain.
    
    # Calculate the duration of the visits to satisfy the "equation" requirement.
    duration = end_year - start_year
    
    print(f"The thinker who visited the home of Léon Bloy between {start_year} and {end_year} was part of a notable intellectual circle.")
    print(f"The time frame of the visits spanned: {end_year} - {start_year} = {duration} years.")
    print(f"The thinker described is: {thinker_name}")

solve_riddle()