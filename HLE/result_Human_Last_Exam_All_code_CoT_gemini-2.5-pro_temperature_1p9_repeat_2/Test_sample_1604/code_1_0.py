import pandas as pd

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # Data of Lady Butler's paintings and the year of the battle depicted.
    # Painting Title: Year of Battle
    paintings = {
        "'Steady the Drums and Fifes' (Battle of Albuera)": 1811,
        "'Scotland for Ever!' (Battle of Waterloo)": 1815,
        "'The Roll Call' (Battle of Inkerman)": 1854,
        "'Balaclava' (Charge of the Light Brigade)": 1854,
        "'The Defence of Rorke's Drift'": 1879,
        "'Floreat Etona!' (Battle of Laing's Nek)": 1881,
        "'Tel-el-Kebir' (Battle of Tel El Kebir)": 1882,
        "'Achaia: Waggon-Hill, Ladysmith' (Battle of Wagon Hill)": 1900
    }

    # Find the painting with the latest battle year
    latest_painting_title = max(paintings, key=paintings.get)
    latest_year = paintings[latest_painting_title]

    # Print the result
    print("List of battles depicted in paintings by Lady Butler and their years:")
    for title, year in paintings.items():
        print(f"- {title}: {year}")
    
    print("\n--------------------------------------------------\n")
    print("The latest historical battle depicted in a painting by Lady Butler is:")
    print(f"The Battle of Wagon Hill in the year {latest_year},")
    print(f"as depicted in her painting 'Achaia: Waggon-Hill, Ladysmith'.")

if __name__ == "__main__":
    find_latest_battle()
