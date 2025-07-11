def solve_historical_riddle():
    """
    This function solves the riddle based on the provided clues.
    """
    # The clues point to a specific historical figure celebrated by Andrei Voznesensky.
    # 1. The poet Voznesensky is famous for the libretto of the rock opera "Juno and Avos".
    # 2. The opera's protagonist is a Russian explorer, fitting the "sailor" description.
    # 3. This man died in Krasnoyarsk, and his grave was lost for many years, inspiring
    #    the line about a missing grave.
    # 4. Efforts to find and re-dedicate the grave began in the late 1980s.
    
    # This person is Nikolai Rezanov.
    last_name = "Rezanov"
    
    explanation = f"""
The person in question is Nikolai Rezanov.

The clues connect to him in the following ways:
- Andrei Voznesensky wrote the poetry for the famous rock opera 'Juno and Avos', which tells Rezanov's story.
- Rezanov was a Russian diplomat and explorer who sailed to Alaska and California, matching the "joy-discovering sailor" description.
- He died in Krasnoyarsk in 1807. His grave was destroyed after the Russian Revolution, and its location was lost for decades, which is why Voznesensky's poem asks about the location of his grave.
- The site of his grave was rediscovered through historical investigation in the late 1980s.

The last name is:
    """
    
    print(explanation)
    print(last_name)

solve_historical_riddle()