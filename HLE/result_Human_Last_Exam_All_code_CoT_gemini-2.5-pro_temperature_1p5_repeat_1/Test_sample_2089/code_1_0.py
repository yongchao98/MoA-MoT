# The film in question is "Collateral".
# Data for the film:
# Release Year: 2004
# Oscar Nominations: 2 (Best Supporting Actor, Best Film Editing)

def calculate_impact_score():
    """
    Calculates a fictional score based on a film's release year and Oscar nominations.
    """
    release_year = 2004
    base_year = 2000
    nominations = 2

    # Calculate the score
    impact_score = (release_year - base_year) * nominations

    # Print the details and the full equation as requested
    print("Film: Collateral")
    print(f"Release Year: {release_year}")
    print(f"Oscar Nominations: {nominations}")
    print("\nCalculating a fictional 'Cinematic Impact Score'...")
    # The final equation with each number printed out
    print(f"Equation: ({release_year} - {base_year}) * {nominations} = {impact_score}")

# Run the calculation
calculate_impact_score()