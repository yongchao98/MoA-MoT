def count_mlb_triple_crowns():
    """
    Calculates and prints the total number of Triple Crown winners in MLB history.
    This includes both batting and pitching Triple Crowns.
    """
    # Historical number of times the batting Triple Crown has been won.
    batting_triple_crown_seasons = 18

    # Historical number of times the pitching Triple Crown has been won.
    pitching_triple_crown_seasons = 36

    # Calculate the total.
    total_seasons = batting_triple_crown_seasons + pitching_triple_crown_seasons

    print("MLB Triple Crown Winners Calculation:")
    print(f"Number of Batting Triple Crown seasons: {batting_triple_crown_seasons}")
    print(f"Number of Pitching Triple Crown seasons: {pitching_triple_crown_seasons}")
    print("-" * 30)
    print("Total Triple Crown seasons in MLB history:")
    # The final equation with each number printed
    print(f"{batting_triple_crown_seasons} (batting) + {pitching_triple_crown_seasons} (pitching) = {total_seasons}")

count_mlb_triple_crowns()