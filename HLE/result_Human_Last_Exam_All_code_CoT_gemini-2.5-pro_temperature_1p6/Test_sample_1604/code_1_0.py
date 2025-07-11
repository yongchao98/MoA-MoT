import operator

def solve():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the years of the battles they depict.
    # Note: Some paintings depict events from longer wars; the year of the specific event or end of the war is used.
    battles_by_butler = {
        "'The 28th Regiment at Quatre Bras' (depicting the Battle of Quatre Bras)": 1815,
        "'Scotland for Ever!' (depicting the Battle of Waterloo)": 1815,
        "'The Roll Call' (depicting the aftermath of the Battle of Inkerman)": 1854,
        "'Balaclava' (depicting the aftermath of the Battle of Balaclava)": 1854,
        "'The Defence of Rorke\'s Drift' (depicting the Battle of Rorke's Drift)": 1879,
        "'The Camel Corps' (depicting a scene from the Gordon Relief Expedition)": 1885,
        "'Within Sound of the Guns' (depicting a scene from the Second Boer War)": 1902
    }

    # Find the painting with the latest date
    latest_painting = max(battles_by_butler.items(), key=operator.itemgetter(1))
    
    # Get all the years and sort them in descending order for the equation
    all_years = sorted(list(set(battles_by_butler.values())), reverse=True)
    
    # Print the result
    print("The latest historical event depicted in a known painting by Lady Butler is from the year {}.".format(latest_painting[1]))
    print("The painting is {}.".format(latest_painting[0]))
    print("\nComparison of the years of the events depicted:")
    
    # Print the equation showing the comparison of the years
    equation = " > ".join(map(str, all_years))
    print(equation)

solve()
<<<Within Sound of the Guns>>>