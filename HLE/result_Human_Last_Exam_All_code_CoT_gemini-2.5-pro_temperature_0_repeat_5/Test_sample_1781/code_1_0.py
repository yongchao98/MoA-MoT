def find_best_historical_answer():
    """
    This script analyzes a historical question about the French monarchy
    and evaluates the provided multiple-choice answers to find the best fit.
    """

    # Historical facts identified from the question
    monarch = "Philip II Augustus"
    event = "Shift from 'King of the Franks' to 'King of France'"
    year_of_event_start = 1190
    year_of_monarchs_death = 1223
    epithet = "Augustus"
    source_of_epithet = "Rigord"

    # Print the analysis step-by-step
    print(f"The question concerns the French monarch {monarch}.")
    print(f"The key event, the shift to a territorial monarchy, began around the year {year_of_event_start}.")
    print(f"The monarch's epithet, '{epithet}', was provided by his biographer, {source_of_epithet}.")
    print(f"The monarch's reign, which contained this transformation, ended with his death in {year_of_monarchs_death}.")
    print("\nEvaluating the choices:")
    print(f"A. 1190, Suetonius: Incorrect person.")
    print(f"B. 1250, Joinville: Incorrect year and person.")
    print(f"C. {year_of_monarchs_death}, {source_of_epithet}: This choice correctly identifies the biographer, Rigord. The year is the end of the monarch's reign, making it the most plausible option.")
    print(f"D. 1789, Voltaire: Incorrect year and person.")
    print(f"E. 1190, Baldwin: Incorrect person (a modern historian, not the source).")
    print("\nConclusion: Choice C is the best fit.")

find_best_historical_answer()