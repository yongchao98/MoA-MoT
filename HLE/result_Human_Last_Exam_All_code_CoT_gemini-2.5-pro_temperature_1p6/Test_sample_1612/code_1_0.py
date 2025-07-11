def count_mlb_triple_crown_winners():
    """
    This function lists and counts the unique players who have won the 
    MLB batting Triple Crown.
    """
    
    # A list of unique players who have won the MLB batting Triple Crown.
    # Note: Rogers Hornsby and Ted Williams each won twice, but are counted once as unique winners.
    winners = [
        # National League Winners
        "Paul Hines",
        "Tip O'Neill",
        "Hugh Duffy",
        "Heinie Zimmerman",
        "Rogers Hornsby",
        "Chuck Klein",
        "Joe Medwick",
        # American League Winners
        "Nap Lajoie",
        "Ty Cobb",
        "Jimmie Foxx",
        "Lou Gehrig",
        "Ted Williams",
        "Mickey Mantle",
        "Frank Robinson",
        "Carl Yastrzemski",
        "Miguel Cabrera"
    ]
    
    total_winners = len(winners)
    
    print("The unique players who have won the MLB Triple Crown are:")
    
    # We will list each player, as each one is a component of the final count.
    for i, name in enumerate(winners):
        print(f"{i + 1}. {name}")
        
    print(f"\nTotal number of unique Triple Crown winners: {total_winners}")

if __name__ == "__main__":
    count_mlb_triple_crown_winners()