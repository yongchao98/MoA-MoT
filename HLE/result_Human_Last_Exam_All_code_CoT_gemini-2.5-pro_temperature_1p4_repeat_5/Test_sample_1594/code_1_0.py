def find_building_shape():
    """
    This function explains the process of identifying the letter formed by
    the buildings of two highschools in Pisa.
    """
    school_1 = "Istituto Galilei-Pacinotti"
    school_2 = "Liceo Scientifico Ulisse Dini"
    city = "Pisa, Italy"

    print(f"Step 1: Locate the schools '{school_1}' and '{school_2}' in {city}.")
    print("Action: Using a mapping service, we find they are adjacent on Via Benedetto Croce.")
    print("-" * 20)

    print("Step 2: Analyze the combined shape of the buildings from a satellite view.")
    print("Observation 1: The building of 'Istituto Galilei-Pacinotti' forms a long, straight, vertical line.")
    print("Observation 2: The building of 'Liceo Scientifico Ulisse Dini' is C-shaped and attached to the first building.")
    print("-" * 20)

    print("Step 3: Combine the shapes to identify the letter.")
    print("Conclusion: The long building acts as the spine, and the C-shaped building provides the three horizontal arms.")
    print("The resulting letter is 'E'.")

if __name__ == "__main__":
    find_building_shape()
