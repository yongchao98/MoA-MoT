def find_vessel_name():
    """
    This function stores the information about the Puntland-Hart Security
    coastguard initiative and prints the name of the centerpiece vessel.
    """
    # Information from the user's query
    state = "Puntland"
    company = "Hart Security Maritime Services Ltd"
    year = 1999

    # The result of the historical lookup
    vessel_name = "MV Pontresina"

    # Print the answer in a full sentence
    print(f"The centerpiece of the coastguard fleet formed by {state} and {company} in {year} was the {vessel_name}.")

if __name__ == "__main__":
    find_vessel_name()