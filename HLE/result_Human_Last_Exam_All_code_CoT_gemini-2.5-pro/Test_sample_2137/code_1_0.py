def solve_riddle():
    """
    This script breaks down the clues to find the correct project.
    """
    
    # Clue 1: The songs and their author.
    songs = ['Venus in Furs', 'Sister Ray', 'Lady Godiva\'s Operation']
    band = "The Velvet Underground"
    principal_songwriter = "Lou Reed"
    
    print("Step 1: Analyzing the songs.")
    print(f"The songs {songs} are by the band '{band}'.")
    print(f"The principal songwriter is {principal_songwriter}.")
    print("-" * 30)

    # Clue 2: Thematic connections and collaborators.
    themes_associated_with = "Andy Warhol and his studio, The Factory"
    key_collaborator = "John Cale"

    print("Step 2: Identifying the central theme and collaborator.")
    print(f"The themes of these songs are deeply connected to {themes_associated_with}, who was the band's manager and producer.")
    print(f"A major musical partner of {principal_songwriter} from this era was {key_collaborator}, co-founder of the band.")
    print("-" * 30)
    
    # Step 3: Synthesizing the clues to find the project.
    print("Step 3: Finding the project.")
    print(f"We are looking for a musical project created as a collaboration between {principal_songwriter} and {key_collaborator} that is about {themes_associated_with.split(' and')[0]}.")
    
    # "Drella" was a nickname for Andy Warhol.
    project_nickname = "Drella"
    final_project = "Songs for Drella"
    
    print(f"The 1990 album '{final_project}' is a song cycle by Reed and Cale dedicated to the life of Andy Warhol, nicknamed '{project_nickname}'.")
    print("-" * 30)
    
    # Final conclusion
    print(f"The correct project is: {final_project}")

solve_riddle()