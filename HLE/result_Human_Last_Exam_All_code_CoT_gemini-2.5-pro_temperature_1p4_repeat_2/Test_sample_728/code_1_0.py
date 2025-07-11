def evaluate_bechdel_test():
    """
    Evaluates a list of films against the Bechdel Test and prints the results.
    """
    films = [
        {'id': 'a', 'title': 'Igla', 'year': 1988, 'rule1': True, 'rule2': True, 'rule3': False, 'reason': "Fails rule 3: Dina and Spartak only talk about a man (Moro)."},
        {'id': 'b', 'title': 'Ghost Town', 'year': 2008, 'rule1': True, 'rule2': False, 'rule3': False, 'reason': "Fails rule 2: The two named women (Gwen and Marjorie) never speak to each other."},
        {'id': 'c', 'title': 'Girls Will Be Girls', 'year': 2003, 'rule1': True, 'rule2': True, 'rule3': True, 'reason': "Passes: Evie and Coco have conversations about their careers and lives."},
        {'id': 'd', 'title': 'War Dogs', 'year': 2016, 'rule1': True, 'rule2': False, 'rule3': False, 'reason': "Fails rule 2: The two named women (Iz and Vanessa) do not speak to each other."},
        {'id': 'e', 'title': 'Slither', 'year': 2006, 'rule1': True, 'rule2': True, 'rule3': False, 'reason': "Fails rule 3: Starla and Kylie's conversation is about a man."},
        {'id': 'f', 'title': 'John Dies at the End', 'year': 2012, 'rule1': True, 'rule2': True, 'rule3': False, 'reason': "Fails rule 3: The conversation between Amy and Jennifer is about a man (John)."},
        {'id': 'g', 'title': 'Man Who Knew Too Much, The', 'year': 1934, 'rule1': True, 'rule2': True, 'rule3': False, 'reason': "Fails rule 3: Jill and Agnes talk to each other, but only about a man."},
        {'id': 'h', 'title': 'Ladies In Retirement', 'year': 1941, 'rule1': True, 'rule2': True, 'rule3': True, 'reason': "Passes: Multiple female characters discuss topics like the house, money, and their situation."},
        {'id': 'i', 'title': 'Internship, The', 'year': 2013, 'rule1': True, 'rule2': True, 'rule3': False, 'reason': "Fails rule 3: Dana and Marielena's conversation is about a man."},
        {'id': 'j', 'title': 'Tinker Bell and the Lost Treasure', 'year': 2009, 'rule1': True, 'rule2': True, 'rule3': True, 'reason': "Passes: The female characters discuss the moonstone, their fairy duties, and the scepter."},
    ]

    passing_films_ids = []

    print("Bechdel Test Evaluation:")
    print("------------------------")
    for film in films:
        passes_test = film['rule1'] and film['rule2'] and film['rule3']
        status = "Passes" if passes_test else "Fails"
        print(f"{film['id']}) {film['title']} ({film['year']}): {status}")
        print(f"   Reason: {film['reason']}\n")
        if passes_test:
            passing_films_ids.append(film['id'])

    # Format the final answer as a comma-separated string
    final_answer = ",".join(passing_films_ids)
    print("Films that pass all three conditions of the Bechdel Test:")
    print(final_answer)

evaluate_bechdel_test()
<<<c,h,j>>>