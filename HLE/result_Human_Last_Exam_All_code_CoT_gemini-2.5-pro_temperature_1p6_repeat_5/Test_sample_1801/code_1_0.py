def find_laureates():
    """
    This function identifies and prints the names of the two Nobel laureates
    featured in a 2000 Enron commercial.
    """
    laureate1_name = "Paul Krugman"
    laureate1_year = 2008

    laureate2_name = "William F. Sharpe"
    laureate2_year = 1990

    print("The two Nobel laureates featured in an infamous 2000 commercial for Enron were:")
    print(f"1. {laureate1_name} (Awarded the Nobel Memorial Prize in Economic Sciences in {laureate1_year})")
    print(f"2. {laureate2_name} (Awarded the Nobel Memorial Prize in Economic Sciences in {laureate2_year})")

if __name__ == "__main__":
    find_laureates()