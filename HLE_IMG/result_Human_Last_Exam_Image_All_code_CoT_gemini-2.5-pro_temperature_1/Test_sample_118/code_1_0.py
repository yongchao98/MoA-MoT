def solve_zurich_path():
    """
    This script determines the optimal path from Billoweg to Laubiweg
    based on the provided criteria by comparing the best 0-exchange path
    with the best 1-exchange path.
    """

    # --- Path Data Derived from Map Analysis ---

    # Path 1: The direct route with 0 exchanges on Line 72.
    path1 = {
        "name": "Direct on Line 72",
        "exchanges": 0,
        "stations": 22, # Manually counted from the map
        "lines": [72],
        "exchange_stations": []
    }

    # Path 2: The best identified alternative with 1 exchange.
    # This path involves taking Line 72 and switching to Line 33 at Haldenegg.
    path2 = {
        "name": "Exchange at Haldenegg",
        "exchanges": 1,
        "stations": 21, # 15 stops on Line 72 + 7 stops on Line 33 - 1 for the exchange station
        "lines": [72, 33],
        "exchange_stations": ["Haldenegg"]
    }

    # --- Decision Logic Based on Criteria ---

    # According to the rules, we first compare paths by the number of exchanges (Criterion A).
    if path1["exchanges"] < path2["exchanges"]:
        best_path = path1
    elif path2["exchanges"] < path1["exchanges"]:
        best_path = path2
    else:
        # If exchanges were equal, we would check Criterion B (fewer stations)
        # and then Criterion C (switch lateness). In this case, it's not necessary.
        # For simplicity, we assume path1 is the default in a tie here.
        best_path = path1

    # --- Format and Print the Final Answer ---

    # Construct the final output string based on the best path found.
    # The format is "S - <line(s)> - E; <station_count>".
    if best_path["exchanges"] == 0:
        line_number = best_path["lines"][0]
        station_count = best_path["stations"]
        final_answer = f"S - {line_number} - E; {station_count}"
    else:
        # This formatting is for paths with exchanges, not needed for the optimal path here.
        path_str_parts = ["S"]
        for i, line in enumerate(best_path["lines"]):
            path_str_parts.append(str(line))
            if i < len(best_path["exchange_stations"]):
                path_str_parts.append(best_path["exchange_stations"][i])
        path_str_parts.append("E")
        station_count = best_path["stations"]
        final_answer = f"{' - '.join(path_str_parts)}; {station_count}"

    print("The optimal path is determined by the following logic:")
    print(f"1. A direct path exists on Line {path1['lines'][0]}, with {path1['exchanges']} exchanges and {path1['stations']} stations.")
    print(f"2. Any alternative path requires at least 1 exchange.")
    print("3. Based on Criterion A (preferring fewer exchanges), the direct path is optimal.")
    print("\nFinal Answer:")
    print(final_answer)

solve_zurich_path()